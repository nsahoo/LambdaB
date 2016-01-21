

echo -e "---- JOB START ----"

cmsRun lambdab_MC.py
echo -e "lambdaB mc Ntuple produced."

mv LambdaB.root LambdaB_mc.root
echo -e "renamed the root file."

cmsRun lambdab_Run2012.py
echo -e "lambdaB data Ntuple produced."

mv LambdaB.root LambdaB_data.root
echo -e "renamed the root file."

echo -e "----  JOB FINISH ----"
