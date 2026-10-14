"""Regression checks for finest-solution acceptance and coarse-level evidence."""
import copy
import json
from pathlib import Path
import validate as v


def main():
    recorded=json.loads(Path(__file__).with_name('results.json').read_text())
    expected=recorded['validation_pass']
    v.evaluate(recorded)
    assert recorded['validation_pass']==expected
    # Controlled gate fixture; this never modifies the saved numerical evidence.
    evidence=copy.deepcopy(recorded)
    for row,balance in zip(evidence['meshes'],(.001,.0015,.001)):
        row['bem']['automatic']['power_balance_rel']=balance
    v.evaluate(evidence)
    assert evidence['validation_pass']
    coarse=copy.deepcopy(evidence)
    coarse['meshes'][0]['bem']['automatic']['power_balance_rel']=.03
    v.evaluate(coarse)
    assert not coarse['validation_pass'] and not coarse['gates']['all_level_bem_balance']
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
    reject(lambda x:x['meshes'][1]['bem']['automatic'].update(power_balance_rel=1.),
           'all_level_bem_balance')
    for key,gate in [('linear_residual_rel','bem_residual'),('faraday_residual_rel','faraday'),
                     ('unit_jump_error','unit_jump')]:
        reject(lambda x,k=key:x['meshes'][0]['bem']['automatic'].update({k:1.}), 'mesh_0.'+gate)
    reject(lambda x:x['meshes'][-1]['bem']['carriers'][0]['l2'].update(power_balance_rel=.03),
           'carrier_solution_checks')
    reject(lambda x:x['meshes'][-1]['bem']['widths']['l2'].update(power_width_rel=.03),'carrier_width')
    planning=copy.deepcopy(evidence)
    planning['runs'][0]['elapsed_s']=1e9
    planning['elapsed_s']=1e10
    v.evaluate(planning)
    assert planning['validation_pass']
    assert not planning['observations']['balance_strictly_decreasing']
    assert not {'runtime','total_runtime','balance_monotonic'} & planning['gates'].keys()
    reject(lambda x:x['meshes'].pop(1),'refinement')
    assert v.LIMITS['bem_power_balance']==v.LIMITS['bem_fem_power']==.02
    print(f'PASS: recorded verdict consistency, controlled acceptance and {checks} rejection checks; fixed 2% thresholds')

if __name__=='__main__':
    main()
