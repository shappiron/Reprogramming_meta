# the function tests if the all ultrasets have consistent lengths between sample and time lists
def ultraset_test(us):
    gse = us.keys()
    for g in gse:
        exps = us[g].keys()
        for e in exps:
            assert len(us[g][e][0]) == len(us[g][e][1]), 'Incosistent samples - time lengths'
    print('Everything is Ok!')   