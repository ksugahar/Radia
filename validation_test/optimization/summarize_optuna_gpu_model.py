"""Summarize recovered TPE experiments; raw operational inputs remain private."""
from pathlib import Path
import argparse
import json,statistics,hashlib
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input-directory",required=True,type=Path)
parser.add_argument("--output",required=True,type=Path)
args=parser.parse_args()
data=args.input_directory; p=data
initial=json.loads((data/'run-metadata.json').read_text());fair=json.loads((data/'fair-metadata.json').read_text());assert len(initial)==6 and len(fair)==10
assert all(not m['preflight']['jobs'] and not m['postflight']['jobs'] and not m['interference'] for m in initial+fair)
independent=json.loads((data/'independent-metadata.json').read_text());assert len(independent)==4
assert all(not m['interference'] and not m['preflight']['jobs'] and not m['postflight']['jobs'] for m in independent)
cub=json.loads((data/'cub-metadata.json').read_text());assert len(cub)==4
assert all(not m['interference'] and not m['preflight']['jobs'] and not m['postflight']['jobs'] for m in cub)
configs=[];references={};compared=0
counter_names=['density_launches','selection_launches','cpu_density_calls','cpu_selection_calls','cpu_parallel_density_calls','host_to_device_bytes','device_to_host_bytes','device_allocation_calls','workspace_reuses','avoided_upload_bytes','model_build_launches','model_reuses','model_download_bytes','preparation_launches','cpu_preparation_calls','cpu_cached_normalizer_values','model_radix_sorts']

def count_values(value):return sum(count_values(v) for v in value) if isinstance(value,list) else 1
def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()
def median(values):return statistics.median(values)
settings=[('initial','cpu',16,'baseline'),('initial','gpu',16,'baseline'),('initial','gpu',16,'candidate')]+[('fair','cpu',t,'candidate') for t in [1,4,8,16]]+[('fair','gpu',16,'candidate')]+[('cub',mode,16,'candidate') for mode in ['cpu','gpu']]
for stage,mode,threads,variant in settings:
 config={'stage':stage,'mode':mode,'threads':threads,'variant':variant,'phases':{}}
 receipt={'initial':initial,'fair':fair,'cub':cub}[stage]
 config['engine_startup_s']=[m['startup_s'] for m in receipt if m['mode']==mode and m['variant']==variant and m['threads']==threads]
 for phase in ['seq','batch','reuse','scaling']:
  if (stage=='fair' and threads!=16 or stage=='cub') and phase in ['batch','reuse']:continue
  pattern=f'{phase}-{mode}-*-{variant}.json' if stage=='initial' else f'{phase}-{stage}-{mode}-{threads}-*-{variant}.json'
  records=[json.loads(f.read_text()) for f in sorted(data.glob(pattern))];assert len(records)==2,(pattern,len(records))
  cases=[]
  for index,ref in enumerate(records[0]['entries']):
   key=(phase,index)
   if key not in references:references[key]=ref
   reference=references[key]
   for r in records:
    entry=r['entries'][index]
    fields=['proposals'] if phase=='scaling' else ['proposals','values','best_curve']
    for field in fields:assert entry[field]==reference[field],(stage,mode,threads,variant,phase,index,field)
    compared+=count_values(entry['proposals'])
   entries=[r['entries'][index] for r in records]
   if phase=='scaling':
    samples=[sum(row)*1000 for e in entries for row in e['samples_s'][4:]]
    cases.append({'history':ref['history'],'dimensions':ref['dimensions'],'warmed_median_ms':median(samples),'warmed_samples_ms':samples,'proposal_sha256':digest(ref['proposals']),'first_trial_s':[e['first_trial_s'] for e in entries],'memory_before_bytes':[e['matlab_memory_before_probes_bytes'] for e in entries],'memory_after_bytes':[e['matlab_memory_after_probes_bytes'] for e in entries]})
   else:
    case={'problem':ref['problem'],'seed':ref['seed'],'dimensions':ref['dimensions'],'proposal_sha256':digest(ref['proposals']),'best_curve_sha256':digest(ref['best_curve']),'best_value':ref['best_value'],'best_at_trials':{str(i):ref['best_curve'][i-1] for i in [10,100,256,512,800]},'inference_total_s':[e['inference_seconds'] for e in entries],'objective_total_s':[e['objective_seconds'] for e in entries],'tell_total_s':[e['history_update_seconds'] for e in entries],'first_proposal_s':[e['first_proposal_seconds'] for e in entries],'trial11_inference_s':[sum(e['samples_s'][10][:2]) for e in entries],'memory_before_bytes':[e['matlab_memory_before_bytes'] for e in entries],'memory_after_bytes':[e['matlab_memory_after_bytes'] for e in entries]}
    case['native_counter_deltas']=[{name:e['native_after'].get(name,0)-e['native_before'].get(name,0) for name in counter_names} for e in entries]
    case['native_workspace_bytes']=[e['native_after'].get('cached_device_bytes',0) for e in entries]
    case['native_peak_sort_temporary_bytes']=[e['native_after'].get('peak_sort_temporary_bytes',0) for e in entries]
    case['native_peak_arena_bytes']=[e['native_after'].get('peak_temporary_bytes',0) for e in entries]
    case['windows']=[]
    for start,end in [(11,100),(101,256),(257,512),(513,800)]:
     samples=[row for e in entries for row in e['samples_s'][start-1:end]];times=sorted(sum(row[:2])*1000 for row in samples)
     case['windows'].append({'trials':[start,end],'inference_median_ms':median(times),'inference_p10_ms':times[int(.1*len(times))],'inference_p90_ms':times[int(.9*len(times))],'ask_median_ms':median([row[0]*1000 for row in samples]),'suggestions_median_ms':median([row[1]*1000 for row in samples])})
    cases.append(case)
  config['phases'][phase]={'module_configure_s':[r['module_configure_s'] for r in records],'cases':cases}
  if phase=='scaling':
   config['phases'][phase]['native_statistics']=[r['api'] for r in records]
  for r in records:
   apis=[r['api']] if phase=='scaling' else [e['native_after'] for e in r['entries']]
   for a in apis:
    if mode=='cpu':assert a['density_launches']==0 and a['host_to_device_bytes']==0 and a['cuda_device']==-1
    else:assert a['density_launches']>0 and a['cpu_density_calls']==0
    if stage in ['fair','cub'] and mode=='gpu':assert a['model_build_launches']>0 and a['cpu_preparation_calls']==0
    if stage in ['fair','cub'] and mode=='cpu':assert a['cpu_cached_normalizer_values']>0 and a['model_build_launches']==0
 configs.append(config)
summary=[]
for phase in ['seq','batch','reuse']:
 for problem in ['branin','correlated_valley']:
  rows=[]
  for c in configs:
   if phase not in c['phases']:continue
   stem=f'{phase}-{c["mode"]}-*-{c["variant"]}.json' if c['stage']=='initial' else f'{phase}-{c["stage"]}-{c["mode"]}-{c["threads"]}-*-{c["variant"]}.json'
   entries=[e for f in data.glob(stem) for e in json.loads(f.read_text())['entries'] if e['problem']==problem]
   samples=[row for e in entries for row in e['samples_s'][512:800]]
   rows.append({'stage':c['stage'],'mode':c['mode'],'threads':c['threads'],'variant':c['variant'],'inference_ms':median([sum(r[:2])*1000 for r in samples]),'ask_ms':median([r[0]*1000 for r in samples]),'suggestions_ms':median([r[1]*1000 for r in samples])})
  final_stage='cub' if phase=='seq' else 'fair'
  cpu=min((r for r in rows if r['mode']=='cpu' and (r['stage']==final_stage or phase=='seq' and r['stage']=='fair' and r['threads']!=16)),key=lambda r:r['inference_ms']);gpu=next(r for r in rows if r['stage']==final_stage and r['mode']=='gpu');old=next(r for r in rows if r['stage']=='initial' and r['mode']=='gpu' and r['variant']=='baseline')
  summary.append({'phase':phase,'gpu_stage':final_stage,'cpu_stage':cpu['stage'],'problem':problem,'warmed_trials':[513,800],'best_measured_cpu_threads':cpu['threads'],'cpu_inference_ms':cpu['inference_ms'],'gpu_inference_ms':gpu['inference_ms'],'previous_gpu_inference_ms':old['inference_ms'],'gpu_speedup_over_best_measured_cpu':cpu['inference_ms']/gpu['inference_ms'],'gpu_speedup_over_previous_gpu':old['inference_ms']/gpu['inference_ms'],'cpu_ask_ms':cpu['ask_ms'],'gpu_ask_ms':gpu['ask_ms'],'cpu_suggestions_ms':cpu['suggestions_ms'],'gpu_suggestions_ms':gpu['suggestions_ms']})
scaling=[]
for i,history in enumerate([24,1000,5000]):
 cpus=[c for c in configs if c['mode']=='cpu' and (c['stage']=='cub' or c['stage']=='fair' and c['threads']!=16)];cpu=min(cpus,key=lambda c:c['phases']['scaling']['cases'][i]['warmed_median_ms']);gpu=next(c for c in configs if c['stage']=='cub' and c['mode']=='gpu');old=next(c for c in configs if c['stage']=='initial' and c['mode']=='gpu' and c['variant']=='baseline');times=[c['phases']['scaling']['cases'][i]['warmed_median_ms'] for c in [cpu,gpu,old]]
 scaling.append({'history':history,'dimensions':40,'cpu_stage':cpu['stage'],'gpu_stage':'cub','best_measured_cpu_threads':cpu['threads'],'cpu_ms':times[0],'gpu_ms':times[1],'previous_gpu_ms':times[2],'gpu_speedup_over_cpu':times[0]/times[1],'gpu_speedup_over_previous_gpu':times[2]/times[1]})
independent_configs=[]
for stage,mode in [('previous_fair','cpu'),('previous_fair','gpu'),('cub','cpu'),('cub','gpu')]:
 rs=[json.loads(f.read_text()) for f in sorted(data.glob(f'ind-cub-{mode}-16-*-candidate.json' if stage=='cub' else f'independent-{mode}-16-*-candidate.json'))];assert len(rs)==2
 entries=[]
 for i,ref in enumerate(rs[0]['entries']):
  key=('independent',i)
  if key not in references:references[key]=ref
  for r in rs:
   assert r['entries'][i]['proposals']==references[key]['proposals']
   compared+=count_values(ref['proposals'])
  samples=[sum(row)*1000 for r in rs for row in r['entries'][i]['samples_s'][4:]]
  entries.append({'history':ref['history'],'dimensions':ref['dimensions'],'warmed_median_ms':median(samples),'warmed_samples_ms':samples,'first_trial_s':[r['entries'][i]['first_trial_s'] for r in rs],'proposal_sha256':digest(ref['proposals'])})
 independent_configs.append({'mode':mode,'threads':16,'stage':stage,'module_configure_s':[r['module_configure_s'] for r in rs],'entries':entries,'native_statistics':[r['api'] for r in rs]})
independent_comparison=[]
for cpu,gpu in zip(independent_configs[2]['entries'],independent_configs[3]['entries']):
 independent_comparison.append({'history':cpu['history'],'dimensions':cpu['dimensions'],'cpu_ms':cpu['warmed_median_ms'],'gpu_ms':gpu['warmed_median_ms'],'gpu_speedup_over_cpu':cpu['warmed_median_ms']/gpu['warmed_median_ms']})
manifest=json.loads((p/'source-hashes.json').read_text());fairmanifest=json.loads((p/'fair-source-hashes.json').read_text());cubmanifest=json.loads((p/'cub-source-hashes.json').read_text());runtime=json.loads((data/'runtime.json').read_text());runtime['cpu']=initial[0]['preflight']['cpu'];runtime['minimum_free_memory_kib']=min(m['preflight']['memory']['FreePhysicalMemory'] for m in initial+fair+cub)
result={'schema':'radia.validation.optuna-gpu-model-performance.v1','date':'2026-10-13','host_role':'on-demand-a100-timing','formal_performance_acceptance':True,'final_source_commit':cubmanifest['candidate'],'gpu_model_commit':manifest['candidate'],'fair_cpu_commit':fairmanifest['candidate'],'source_archives':{'initial':manifest['archive_sha256'],'fair':fairmanifest['archive_sha256'],'final':cubmanifest['archive_sha256']},'binaries':{'previous':manifest['binary_sha256']['baseline'],'initial_gpu_model':manifest['binary_sha256']['candidate'],'fair_cpu_and_gpu':fairmanifest['binary_sha256'],'final_cpu_and_gpu':cubmanifest['binary_sha256']},'runtime':runtime,'design':{'fresh_engine_sessions':24,'final_order':[[m['mode'],m['threads']] for m in cub],'initial_order':[[m['mode'],m['threads'],m['variant']] for m in initial],'fair_order':[[m['mode'],m['threads'],m['variant']] for m in fair],'candidates':24,'multivariate':True,'seeds':[11,37],'closed_loop_trials':800,'phase_settings':{'seq':{'batch_size':1,'constant_liar':True},'batch':{'batch_size':8,'constant_liar':True},'reuse':{'batch_size':8,'constant_liar':False}},'all_history_retained':True,'storage_enabled':False,'scaling_repeats':13,'scaling_discarded_initial_repeats':4,'scaling_warmed_samples_per_configuration_case':18,'foreign_job_monitor_interval_s':10,'no_foreign_jobs_detected':True,'all_sessions_closed':True,'all_proposals_objectives_and_convergence_curves_exactly_equal':True,'compared_proposal_values':compared},'configurations':configs,'closed_loop_comparison':summary,'fixed_history_comparison':scaling,'independent_configurations':independent_configs,'independent_comparison':independent_comparison,'limitations':['Branin and a 40D correlated quadratic valley are public optimization benchmarks, not engineering CAE acceptance.','ask includes model preparation/sampling and Study setup; suggestions include parameter retrieval/recording. API timings are not kernel-only timings.','CPU 1/4/8/16 are screened for sequential and scaling cases; batch comparisons use CPU/GPU 16 only.','Initial, fair and final CUB/parallel-exponential blocks are separate experiments. Final CPU16/GPU16 use the same binary in ABBA order. CPU1/4/8 screening uses the previous binary; the final changes affect GPU kernels only. Batch/reuse results remain on the previous fair binary and are not final-binary performance claims.','Cold initialization and per-seed totals are recorded; warmed summaries pool trials 513-800.','Only the good model sampling parameters return to the ordered host RNG. Mixture weights, history selection/gathering and domain checks remain host work.','CPU also reuses constructor normalizers and prepares log-sigma once; no repeated-normalizer CPU baseline is used for the best-CPU speedup.','Before/after MATLAB memory, retained device capacity and transfers are measured; whole-process peak memory and isolated transfer seconds are not.','Source-built development artifacts, not a release or a Python/DB, Rustuna or BoTorch speed comparison.']}
args.output.write_text(json.dumps(result,indent=2)+'\n')
for row in summary+scaling+independent_comparison:print(row)
