from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]

def test_paired_costs_and_confidence_intervals():
    df=pd.read_csv(ROOT/'figure4/figure4_fivepoint_SE_CI95.csv')
    old=pd.read_csv(ROOT/'results/controller_performance.csv')
    assert len(df)==73
    assert set(np.round(df.q.unique(),2))=={.92,.96,.98}
    for (_,_,q), group in df.groupby(['a1','a2','q']):
        assert len(np.unique(group.a))==5
        if np.isclose(q,.98):
            st=group[group.policy=='stationary']
            assert len(st)==3
    assert np.all(df.sample_SE>=0)
    assert np.all(df.bootstrap_CI95_low<=df.bootstrap_CI95_high)
    assert np.all(df.sample_SE >= 0)
    raw=np.load(ROOT/'figure4/figure4_trajectory_costs.npz')
    idx=pd.read_csv(ROOT/'figure4/figure4_trajectory_cost_index.csv')
    assert len(idx)==len(df)
    for (i,row) in df.reset_index(drop=True).iterrows():
        values=raw[idx.iloc[i].npz_key]
        assert len(values)==1200
        assert np.isclose(np.mean(values),row.mean_cost,rtol=1e-12,atol=1e-12)
        assert np.isclose(np.std(values,ddof=1)/np.sqrt(len(values)),row.sample_SE,rtol=1e-11,atol=1e-12)
    merged=old.merge(df,left_on=['a1','a2','q','a','controller'],right_on=['a1','a2','q','a','policy'])
    assert len(merged)==43
    np.testing.assert_allclose(merged.mean_J_MC,merged.mean_cost,rtol=5e-12,atol=1e-12)
