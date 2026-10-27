import argparse,contextlib,hashlib,io,json,sys,time,warnings
from pathlib import Path
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--source',type=Path);p.add_argument('--output',type=Path);args=p.parse_args()
sys.path.insert(0,str(args.source));sys.path.insert(0,str(args.source/'csdid/test_csdid'))
from csdid.att_gt import ATTgt
from helpers import build_sim_data

def canonical(value):
 # aggte's call dictionary stores its owning fit object. Encode the backreference
 # by class, while hashing that object's complete MP/results separately.
 if isinstance(value,ATTgt):return {'fit_reference':'csdid.att_gt.ATTgt'}
 if isinstance(value,pd.DataFrame):
  return {'frame':canonical(list(value.columns)),'index':canonical(value.index.to_numpy()),'columns':[canonical(value[c].to_numpy())for c in value.columns],'dtypes':[str(x)for x in value.dtypes]}
 if isinstance(value,pd.Series):return {'series':canonical(value.to_numpy()),'index':canonical(value.index.to_numpy())}
 if isinstance(value,np.ndarray):
  return {'dtype':str(value.dtype),'shape':value.shape,'bytes':value.tobytes().hex()}if value.dtype.kind!='O'else {'dtype':'object','shape':value.shape,'values':canonical(value.tolist())}
 if isinstance(value,dict):return {'dict':[[canonical(k),canonical(v)]for k,v in sorted(value.items(),key=lambda kv:repr(kv[0]))]}
 if isinstance(value,(tuple,list)):return {type(value).__name__:[canonical(x)for x in value]}
 if isinstance(value,(np.floating,float)):return {'float':float(value).hex()}
 if isinstance(value,(np.integer,int))and not isinstance(value,(bool,np.bool_)):return int(value)
 if isinstance(value,(np.bool_,bool)):return bool(value)
 if value is None or isinstance(value,str):return value
 raise TypeError(type(value))

def digest(value):return hashlib.sha256(json.dumps(canonical(value),sort_keys=True,separators=(',',':')).encode()).hexdigest()

base=build_sim_data(n=160,time_periods=5,seed=19371)
base['wt']=np.exp(.2*base.X+.1*np.sin(base.id)+.03*base.period)
base.loc[base.id%17==0,'wt']=0.
records=[];started=time.monotonic()
for scenario in ['panel','varying','rc','unbalanced']:
 for method in ['dr','ipw','reg']:
  for faster in [False,True]:
   for scale in [1e-50,1e-10,1.,1e10,1e50]:
    for bootstrap in [False,True]:
     data=base.copy();data.wt*=scale
     if scenario=='rc':data['id']=np.arange(len(data))
     if scenario=='unbalanced':data=data.loc[~((data.id%11==0)&data.period.eq(3))].copy()
     key=f'{scenario}/{method}/{faster}/{scale}/{bootstrap}'
     np.random.seed(919273)
     with warnings.catch_warnings(record=True)as caught,contextlib.redirect_stdout(io.StringIO()):
      warnings.simplefilter('always')
      f=ATTgt(yname='Y',tname='period',idname='id',gname='G',data=data,panel=scenario!='rc',allow_unbalanced_panel=scenario=='unbalanced',weights_name='wt',fix_weights='varying'if scenario=='varying'else None,xformla='~X',faster_mode=faster,cband=False,biters=199).fit(method,bstrap=bootstrap)
      hashes={'MP':digest(f.MP),'results':digest(f.results)}
      for aggregation in ['simple','dynamic','group','calendar']:
       f.aggte(aggregation,bstrap=bootstrap,cband=False);hashes[aggregation]=digest(f.atte)
     hashes['warnings']=digest([(type(w.message).__name__,str(w.message))for w in caught]);hashes['RNG']=digest(np.random.get_state())
     records.append({'workflow':key,'sha256':hashes})
report={'source':str(args.source),'workflows':len(records),'aggregations':4*len(records),'elapsed_seconds':time.monotonic()-started,'records':records}
args.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items()if k!='records'}))
