import csv, json, sys, statistics as st
survey, emp = sys.argv[1], sys.argv[2]
ids=set(r["sequence"]+"/"+r["charge"] for r in csv.DictReader(open(emp),delimiter="\t"))
rows=[]
for r in csv.reader(open(survey),delimiter="\t"):
    if r[0]=="pair": continue
    rows.append((r[1],float(r[2]),float(r[3]),float(r[4])))
ranked=sorted(range(len(rows)),key=lambda i:-rows[i][3])
rank={rows[i][0]:k for k,i in enumerate(ranked)}
E=[r for r in rows if r[0] in ids]; N=[r for r in rows if r[0] not in ids]
q=lambda v,p: sorted(v)[int(p*(len(v)-1))]
out={"E_c":len(E),
 "E_rank_quantiles":{p:q([rank[r[0]] for r in E],p) for p in (0.1,0.25,0.5,0.75,0.9)},
 "E_target_gt_decoy":sum(r[1]>r[2] for r in E)/len(E),
 "E_target_coelution_median":st.median(r[1] for r in E),"E_decoy_coelution_median":st.median(r[2] for r in E),
 "other_target_coelution_median":st.median(r[1] for r in N),"other_decoy_coelution_median":st.median(r[2] for r in N),
 "S_quantiles_all":{p:q([r[3] for r in rows],p) for p in (0.5,0.75,0.9,0.99)},
 "S_at_rank":{k:rows[ranked[k]][3] for k in (100000,200000,300000,500000,800000)}}
print(json.dumps(out,indent=1))
