"""Guards: pagination, attempts and missing evidence cannot become complete CI claims.

Assumes: synthetic API responses only; no credentials/network or numerical models.
Durations describe recorded intervals, not performance acceptance criteria.
"""
import importlib.util
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('ci_metrics', ROOT/'tools/collect_tests_design_metrics.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class MetricsContracts(unittest.TestCase):
    def test_get_is_explicit_read_only_and_does_not_report_stderr(self):
        with patch.object(m.subprocess, 'run', return_value=subprocess.CompletedProcess([],0,'{}','')) as run:
            self.assertEqual(m.gh_get('repos/example/project/pulls'), {})
            argv=run.call_args.args[0]
            self.assertEqual(argv[:5], ['gh','api','--method','GET','-H'])
        with patch.object(m.subprocess, 'run', return_value=subprocess.CompletedProcess([],1,'','secret')):
            with self.assertRaisesRegex(RuntimeError, 'response not logged') as error:
                m.gh_get('repos/example/project/pulls')
            self.assertNotIn('secret',str(error.exception))

    def test_paginate_more_than_one_hundred_and_hash_every_response(self):
        with TemporaryDirectory() as tmp:
            calls=[]
            def get(url):
                calls.append(url)
                return [{'id':i} for i in range(100)] if url.endswith('&page=1') else [{'id':100}]
            s=m.Snapshot(Path(tmp)/'private',get)
            self.assertEqual(len(s.pages('items')),101)
            self.assertEqual(len(calls),2)
            import hashlib
            for row in s.responses:
                self.assertEqual(hashlib.sha256((s.out/row['path']).read_bytes()).hexdigest(),row['sha256'])

    def test_pagination_does_not_silently_stop_at_search_cap(self):
        from urllib.parse import urlsplit, parse_qs
        with TemporaryDirectory() as tmp:
            def get(url):
                page=int(parse_qs(urlsplit(url).query)['page'][0])
                return {'total_count':1001,'workflow_runs':[
                    {'id':i} for i in range((page-1)*100,min(page*100,1001))]}
            s=m.Snapshot(Path(tmp)/'private',get)
            self.assertEqual(len(s.pages('runs','workflow_runs')),1001)
            self.assertFalse(s.errors)


    def test_partial_page_preserves_prior_observations(self):
        with TemporaryDirectory() as tmp:
            def get(url):
                if url.endswith('&page=1'):return [{'id':i} for i in range(100)]
                raise RuntimeError('secret')
            s=m.Snapshot(Path(tmp)/'private',get)
            self.assertEqual(len(s.pages('items')),100)
            self.assertEqual(len(s.errors),1)
            self.assertNotIn('secret',json.dumps(s.errors))

    def test_repeated_page_and_total_count_mismatch_are_incomplete(self):
        for mode in ('repeat','count'):
            with self.subTest(mode=mode),TemporaryDirectory() as tmp:
                def get(url):
                    return {'jobs':[{'id':i} for i in range(100)],'total_count':200} if mode=='repeat' else {'jobs':[], 'total_count':1}
                s=m.Snapshot(Path(tmp)/'private',get)
                s.pages('jobs','jobs')
                self.assertTrue(s.errors)

    def test_outside_worktree_fresh_output_required(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'.git').mkdir()
            with self.assertRaisesRegex(ValueError,'worktree'):m.Snapshot(root/'output')
        with TemporaryDirectory() as tmp:
            out=Path(tmp)/'private';m.Snapshot(out)
            with self.assertRaises(FileExistsError):m.Snapshot(out)
        with self.assertRaises(ValueError):m.Snapshot(Path('relative'))

    def test_association_never_uses_branch_name(self):
        pr={'number':3,'head':{'sha':'old'},'merge_commit_sha':'merge'}
        self.assertIsNone(m.association(pr,{'head_branch':'same','head_sha':'other'}))
        self.assertEqual(m.association(pr,{'head_sha':'merge'}),'exact_merge_sha')
        self.assertEqual(m.association(pr,{'head_sha':'old'}),'exact_final_head_sha')
        self.assertEqual(m.association(pr,{'pull_requests':[{'number':3}]}),'api_pull_request_link')

    def synthetic(self, url):
        if '/pulls?' in url:
            return [{'id':3,'number':3,'created_at':'2026-01-01T00:00:00Z',
                     'merged_at':'2026-01-02T00:00:00Z','head':{'sha':'head'},'merge_commit_sha':'merge'}]
        if '/actions/runs?' in url:
            return {'total_count':1,'workflow_runs':[{'id':8,'head_sha':'head','run_attempt':2,'event':'pull_request'}]}
        attempt=2 if '/attempts/2' in url else 1
        if '/jobs?' in url:
            return {'total_count':2,'jobs':[
                {'id':100+attempt,'run_id':8,'run_attempt':attempt,'name':'build','status':'completed',
                 'conclusion':'success','started_at':'2026-01-01T00:01:00Z','completed_at':'2026-01-01T00:02:00Z'},
                {'id':200+attempt,'run_id':8,'run_attempt':attempt,'name':'skip','conclusion':'skipped'}]}
        return {'id':8,'run_attempt':attempt,'event':'pull_request','head_sha':'head',
                'created_at':'2026-01-01T00:00:00Z','run_started_at':'2026-01-01T00:01:00Z',
                'updated_at':'2026-01-01T03:00:00Z','status':'completed','conclusion':'success'}

    def test_attempt_jobs_and_elapsed_are_distinct(self):
        with TemporaryDirectory() as tmp:
            s=m.Snapshot(Path(tmp)/'private',self.synthetic)
            result=m.collect(s,'example/project','2026-01-02T00:00:00Z','2026-01-03T00:00:00Z')
            self.assertTrue(result['collection_complete'])
            self.assertEqual([x['attempt'] for x in result['runs']],[1,2])
            self.assertEqual([x['queue_s'] for x in result['runs']],[60,None])
            self.assertTrue(all(x['wall_s'] is None for x in result['runs']))
            self.assertEqual(result['prs'][0]['open_to_merge_s'],86400)
            self.assertEqual(result['prs'][0]['executed_job_count'],2)
            self.assertEqual(result['prs'][0]['skipped_job_count'],2)
            self.assertEqual(result['prs'][0]['observed_sum_job_execution_s'],120)
            self.assertEqual(len(result['jobs']),4)
            schema=json.loads((ROOT/'tools/schemas/tests_design_metrics.schema.json').read_text())
            for key in ('prs','runs','jobs'):
                required=set(schema['properties'][key]['items']['required'])
                for row in result[key]: self.assertEqual(set(row),required)

    def test_mismatched_attempt_jobs_cannot_be_counted(self):
        with TemporaryDirectory() as tmp:
            def get(url):
                data=self.synthetic(url)
                if '/jobs?' in url:
                    for job in data['jobs']:job['run_attempt']=99
                return data
            result=m.collect(m.Snapshot(Path(tmp)/'private',get),'a/b','2026-01-02T00:00:00Z','2026-01-03T00:00:00Z')
            self.assertFalse(result['collection_complete'])
            self.assertEqual(result['jobs'],[])

    def test_no_merged_prs_is_empty_observation_not_zero_latency(self):
        with TemporaryDirectory() as tmp:
            result=m.collect(m.Snapshot(Path(tmp)/'private',lambda url:[]),'a/b','2026-01-02T00:00:00Z','2026-01-03T00:00:00Z')
            self.assertTrue(result['collection_complete'])
            self.assertEqual(result['prs'],[])
            self.assertEqual(result['runs'],[])

    def test_malformed_timestamps_preserve_partial_snapshot(self):
        with TemporaryDirectory() as tmp:
            def get(url):
                data=self.synthetic(url)
                if '/jobs?' in url: data['jobs'][0]['completed_at']='invalid'
                return data
            result=m.collect(m.Snapshot(Path(tmp)/'private',get),'a/b','2026-01-02T00:00:00Z','2026-01-03T00:00:00Z')
            self.assertFalse(result['collection_complete'])
            self.assertIsNone(result['jobs'][0]['duration_s'])
            schema=json.loads((ROOT/'tools/schemas/tests_design_metrics.schema.json').read_text())
            self.assertEqual(set(result),set(schema['required']))


    def test_unknown_and_invalid_time_are_not_zero(self):
        self.assertIsNone(m.elapsed(None,'2026-01-01T00:00:00Z'))
        with self.assertRaises(ValueError):m.elapsed('2026-01-02T00:00:00Z','2026-01-01T00:00:00Z')
        with self.assertRaises(ValueError):m.timestamp('2026-01-01T00:00:00')


if __name__ == '__main__':
    unittest.main()
