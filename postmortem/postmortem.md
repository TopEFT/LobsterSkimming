# Post-mortem reweighting
## Setup
This is based on the [CMGTools repo](https://github.com/sscruz/cmgtools-lite/tree/2aa24e61e5500bfdb490e93902c256be5fc29768?tab=readme-ov-file#set-up-cmssw-and-the-base-git)

Install CMSSW and setup the base git
```
cmsrel CMSSW_10_6_26
cd CMSSW_10_6_26/src
cmsenv
git cms-init
```

Get the central Heppy 80X branch (might not be needed, should check)
```
git remote add cmg-central https://github.com/CERN-PH-CMG/cmg-cmssw.git -f  -t heppy_80X
cp /afs/cern.ch/user/c/cmgtools/public/sparse-checkout_80X_heppy .git/info/sparse-checkout
git checkout -b heppy_80X cmg-central/heppy_80X
```

Get the CMGTools subsystem
```
git clone -o cmg-central https://github.com/CERN-PH-CMG/cmgtools-lite.git -b 80X CMGTools
cd CMGTools
```

Compile CMSSW
```
cd $CMSSW_BASE/src
scram b -j 8
```
