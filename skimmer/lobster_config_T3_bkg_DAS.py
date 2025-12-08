import datetime
import os
import sys
import shutil
import subprocess
import inspect

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, Dataset, ParentDataset, StorageConfiguration, Workflow

sys.path.append(os.path.split(__file__)[0])
from tools.utils import read_cfg

TSTAMP1 = datetime.datetime.now().strftime('%Y%m%d_%H%M')
startingday = datetime.datetime.now().strftime('%y%m%d')

top_dir = subprocess.check_output(["git", "rev-parse", "--show-toplevel"]).decode("utf-8").strip()
print(top_dir)

#sandbox_location = os.path.join(top_dir,"LobsterSkimming/CMSSW_14_0_6")
sandbox_location = "/users/apiccine/work/LobPY3/LobsterSkimming/CMSSW_14_0_6"
#sandbox_location = "/users/apiccine/work/Lobster-With-Conda/LobsterSkimming/CMSSW_14_0_6"

testing = False #True

step = "run3skimsT3"
#ver = "v1"
ver = "v{}".format(startingday)
year = "2023BPix"

tag = f"mc/NAOD_ULv12_lepMVA-run3/{year}/Full"
#cfg_name = f"ND_{year}_background_samples.cfg"
cfg_name = f"central_{year}_ZG_samples.cfg"


cfg_fpath = os.path.join(top_dir,"topeft/input_samples/cfgs",cfg_name)

# Only process json files that match these regexs (empty list matches everything)
match = ['.*ZG_MLL.*\\.json']
#match = ['ZG_MLL-4to50_PTG-10to100\\.json',  'ZG_MLL-50_PTG-600\\.json']
#match = ['ZG_MLL-4to50_PTG-10to100\\.json']
#match = ['.*ZZZ.*\\.json']
#match = ['.*DYJetsToLL_MLL-50.*\\.json'] #, '.*ZZTo4L.*\\.json']
#match = ['.*TWZ_Tto2Q_WtoLNu_Zto2L.*\\.json', '.*DYJetsToLL_MLL-50.*\\.json']
#match = ['.*.json']
#match = ['.*Muon.*22Sep2023\\.json']
# match = ['DoubleEG_F-UL2016\\.json']
# match = ['MuonEG_B-UL2017\\.json']

skim_cut = "nMuon+nElectron+nTau>=2&&Sum\$(Muon_looseId&&Muon_miniPFRelIso_all<0.4&&Muon_sip3d<8)+Sum\$(Electron_miniPFRelIso_all<0.4&&Electron_sip3d<8)+Sum\$(Tau_idDeepTau2018v2p5VSe>1&&Tau_idDeepTau2018v2p5VSmu>0&&Tau_idDeepTau2018v2p5VSjet>1)>=2"

master_label = '{step}_lobPY3_{cfg}_{tstamp}'.format(step=step,cfg=cfg_name.split(".")[0],tstamp=TSTAMP1)
workdir_path = "{path}/{step}/{tag}/{ver}".format(step=step,tag=tag,ver=ver,path="/scratch365/$USER")
plotdir_path = "{path}/{step}/{tag}/{ver}".format(step=step,tag=tag,ver=ver,path="~/www/lobster")
output_path  = "{path}/{step}/{tag}/{ver}".format(step=step,tag=tag,ver=ver,path="/store/user/$USER")

if testing:
    master_label = '{step}_testlobPY3_{tstamp}'.format(step=step,tstamp=TSTAMP1)
    workdir_path = "{path}/{step}/test/lobster_skimtest_{tstamp}".format(step=step,tstamp=TSTAMP1,path="/scratch365/$USER")
    plotdir_path = "{path}/{step}/test/lobster_skimtest_{tstamp}".format(step=step,tstamp=TSTAMP1,path="~/www/lobster")
    output_path  = "{path}/{step}/test/lobster_skimtest_{tstamp}".format(step=step,tstamp=TSTAMP1,path="/store/user/$USER")

# Different xrd src redirectors depending on where the inputs are stored
#xrd_src = "ndcms.crc.nd.edu"            # Use this for accessing samples from the GRID
#xrd_src = "cmsxrootd.fnal.gov"          # Only use this if the ND XCache is giving troubles
xrd_src = "skynet013.crc.nd.edu"         # Use this for accessing samples from ND T3
#xrd_src = "cms-xrd-global.cern.ch"       # global redirector
xrd_dst = "skynet013.crc.nd.edu"

storage_base = StorageConfiguration(
    input=[
        #"file:///cms/cephfs/data/",
        "root://{host}//".format(host=xrd_src)  # Note the extra slash after the hostname
    ],
    output=[
        f"file:///cms/cephfs/data/{output_path}",
        #f"root://{xrd_dst}/{output_path}",
    ],
    #disable_input_streaming=True,
)


storage_cmssw = StorageConfiguration(
    #input=[
    #    #"file:///cms/cephfs/data/",
    #    "root://{host}//".format(host=xrd_src)  # Note the extra slash after the hostname
    #],
    output=[
        f"file:///cms/cephfs/data/{output_path}",
        f"root://{xrd_dst}/{output_path}",
    ],
    disable_input_streaming=False,
)

storage = storage_cmssw
#storage = storage_base

# See tools/utils.py for dict structure of returned object
cfg = read_cfg(cfg_fpath,match=match)
skim_cut = skim_cut.replace(">", "\>").replace("<", "\<").replace(" ", "")

cat = Category(
    name='processing',
    cores=1,
    memory=59000,
    disk=47000,
)

samples = {}

wf = []
for sample in sorted(cfg['jsons']):
    jsn = cfg['jsons'][sample]
    year = jsn['year']
    daspath = jsn['path']
    print(jsn['path'])
    
    if year not in list(samples.keys()):
        samples[year] = []
    samples[year].append(daspath)

    module_name = ''
    if '16' in year:
        module_name = 'lepMVA_2016_preVFP'
    elif '17' in year:
        module_name = 'lepMVA_2017'
    elif '18' in year:
        module_name = 'lepMVA_2018'
    elif '2022' or '2023' in year:
        module_name = 'lepMVA'
    else:
        module_name = 'lepMVA_2016'
            
    ds_cmssw = cmssw.Dataset(
        dataset=jsn['path'],
        lumis_per_task=1,   # Since file_based is true, this will be files_per_task
        file_based=True
    )

    ds = ds_cmssw
    
    cmd = ['python3','skim_wrapper_T3.py']
    cmd.extend(['--cut',skim_cut])
    cmd.extend(['--module',module_name])
    cmd.extend(['--out-dir','.'])
    cmd.extend(['@inputfiles'])

    print("\n\n\n")
    print("Command to execute:", ' '.join(cmd))
    print("\n\n\n")

    skim_wf = Workflow(
        label=sample.replace('-','_'),
        sandbox=cmssw.Sandbox(release=sandbox_location),
        dataset=ds,
        category=cat,
        extra_inputs=['skim_wrapper_T3.py'], #,os.path.join(sandbox_location,'src/PhysicsTools/NanoAODTools/scripts/haddnano.py')],
        outputs=['output.root'],
        command=' '.join(cmd),
        #merge_command='python haddnano.py @outputfiles @inputfiles',
        merge_command='haddnano.py @outputfiles @inputfiles',
        merge_size='537M',
        globaltag=False,    # To avoid autosense crash (can be anything, just not None)
        cleanup_input=False
    )
    wf.extend([skim_wf])
    
config = Config(
    label=master_label,
    workdir=workdir_path,
    plotdir=plotdir_path,
    storage=storage,
    workflows=wf,
    advanced=AdvancedOptions(
	bad_exit_codes=[127, 160],
        log_level=1,
        payload=10,
	osg_version='3.6',
        threshold_for_failure=1000,
	threshold_for_skipping=1000,
    )
)
