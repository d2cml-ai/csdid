# -------------- package
# from .mp import mp_py
# from .pre_process_did import pre_process_did
# from .compute_attgt import compute_att_gt
# from .process_attgt import process_attgt
# -------------- local dev
from mp import mp_py
from pre_process_did import pre_process_did
from compute_attgt import compute_att_gt
from process_attgt import process_attgt

def process_attgt(
    yname,
    tname,
    gname,
    data,
    idname=None,
    xformla=None,
    panel=True,
    allow_unbalanced_panel=False,
    control_group=["nevertrated", "notyettreated"],
    anticipation=0,
    weightsname=None,
    alpha=0.05,
    bstrap=True,
    cband=True,
    cband=True,
    biters=1000,
    clustervars=None,
    est_method="dr",
    base_period="varying",
    print_details=False,
    pl=False,
    cores=1,
  ) -> list:
  dp = pre_process_did(
    yname,
    tname,
    idname,
    gname,
    data,
    panel,
    
  )
