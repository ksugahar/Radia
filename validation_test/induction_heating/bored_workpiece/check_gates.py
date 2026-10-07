"""Regression checks for finest-solution acceptance and coarse-level evidence."""
import copy
import json
from pathlib import Path
import validate as v


def main():
    evidence=json.loads(Path(__file__).with_name('results.json').read_text())
    v.evaluate(evidence)
    assert evidence['validation_pass']
    assert evidence['meshes'][0]['bem']['automatic']['power_balance_rel']>.02
    assert evidence['acceptance_reference']==dict(bem_nodes=6232,fem_ndof=601029)
    checks=0
    def reject(edit,gate):
        nonlocal checks
        x=copy.deepcopy(evidence);edit(x);v.evaluate(x)
        result=x['gates']
        for key in gate.split('.'): result=result[key]
        assert not result and not x['validation_pass'],gate
        checks+=1
    reject(lambda x:x['meshes'][-1]['bem']['automatic'].update(power_balance_rel=.0201),
           'finest_solution.bem_balance')
    reject(lambda x:x['meshes'][0]['fem'].update(P_total=70.),'finest_solution.agreement')
    reject(lambda x:x['meshes'][1]['bem']['automatic'].update(power_balance_rel=.03),
           'balance_monotonic')
    for key,gate in [('linear_residual_rel','bem_residual'),('faraday_residual_rel','faraday'),
                     ('unit_jump_error','unit_jump')]:
        reject(lambda x,k=key:x['meshes'][0]['bem']['automatic'].update({k:1.}), 'mesh_0.'+gate)
    reject(lambda x:x['meshes'][-1]['bem']['carriers'][0]['l2'].update(power_balance_rel=.03),
           'carrier_solution_checks')
    reject(lambda x:x['meshes'][-1]['bem']['widths']['l2'].update(power_width_rel=.03),'carrier_width')
    reject(lambda x:x['runs'][0].update(elapsed_s=1201),'runtime')
    reject(lambda x:x.update(elapsed_s=3601),'total_runtime')
    reject(lambda x:x['meshes'].pop(1),'refinement')
    assert v.LIMITS['bem_power_balance']==v.LIMITS['bem_fem_power']==.02
    print(f'PASS: baseline and {checks} rejection checks; fixed 2% thresholds')

if __name__=='__main__':
    main()
