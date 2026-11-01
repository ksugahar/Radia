"""GPU-only lifetime contracts through real production entries, no fake kernels.

Executed only by the separately admitted native/GPU validation runner. There is
no precision criterion here: holding a second model must not mutate the first,
and shutdown must refuse live owners. Physical reset/address-reuse qualification
remains a separate exclusive-window test, not certified by this contract.
"""
def validate(k):
    d=k.Dimension();d.low=-2;d.high=2;d.step=float('nan');d.minimum=0
    d.choices=0;d.categorical=False;d.log_scale=False
    weights=[0.25]*4
    anchor=k.build_model('gpu',[0.8,0.9,1.0],weights,[d],1.0,False,False)
    first=k.build_model('gpu',[-1.0,-0.9,-0.8],weights,[d],1.0,False,True)
    candidates=[0.0,-1.5,-1.0,-0.5,0.5,1.0,1.5]
    before=k.best_prepared('gpu',candidates,first,anchor)
    second=k.build_model('gpu',[0.8,0.9,1.0],weights,[d],1.0,False,True)
    after=k.best_prepared('gpu',candidates,first,anchor)
    alternative=k.best_prepared('gpu',candidates,second,anchor)
    assert before!=alternative,'Lifecycle inputs do not distinguish the overwritten-model mutation'
    assert after==before,('A later build overwrote a live model',before,after)
    try:k.shutdown()
    except RuntimeError as error:
        assert 'live buffers' in str(error),str(error)
    else:raise AssertionError('Shutdown discarded callbacks with live device owners')
    for model in (first,second,anchor):model.close()
    try:k.best_prepared('gpu',candidates,first,anchor)
    except ValueError as error:assert 'closed' in str(error),str(error)
    else:raise AssertionError('A closed model was accepted')
    return {'pinned_model_unchanged':True,'witness_sensitive':True,
            'shutdown_refused_live_owners':True,'closed_owner_rejected':True}
