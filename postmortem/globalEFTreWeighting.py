from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import writeOutput
import os, sys
import math
import numpy as np
import imp
import tempfile, shutil, stat
import numpy
import json
#from ROOT import TMatrixD, TDecompSVD, TVectorD
#import array

def invert_momenta(p):
    #fortran/C-python do not order table in the same order
    new_p = []
    for i in range(len(p[0])):
        new_p.append([0]*len(p))
    for i, onep in enumerate(p):
        for j, x in enumerate(onep):
            new_p[j][i] = x
    return new_p

def SortPDGs(pdgs):
    return sorted(pdgs[:2]) + sorted(pdgs[2:])


def zboost(part, pboost=[]):
    """Both momenta should be in the same frame.
The boost perform correspond to the boost required to set pboost at
    rest (only z boost applied).
    """
    E = pboost[0]
    pz = pboost[3]
    #beta = pz/E
    gamma = E / math.sqrt(E**2-pz**2)
    gammabeta = pz  / math.sqrt(E**2-pz**2)

    out =  [gamma * part[0] - gammabeta * part[3],
            part[1],
            part[2],
            gamma * part[3] - gammabeta * part[0]]

    if abs(out[3]) < 1e-6 * out[0]:
        out[3] = 0
    return out


def numToString(num):
    return ("%4.2f"%num).replace('.','p').replace('-','m')

class globalEFTreWeighting( Module ):
    def __init__(self, process):

        self.process_dict = {
            'tHq' : "tHq4f_all22WCsStartPtCheckdim6TopMay20GST_run0",
            'tllq' : "tllq4fNoSchanWNoHiggs0p_all22WCsStartPtCheckV2dim6TopMay20GST_run0",
            'ttH' : "ttHJet_all22WCsStartPtCheckdim6TopMay20GST_run0",
            'ttll' : "ttllNuNuJetNoHiggs_all22WCsStartPtCheckdim6TopMay20GST_run0",
            'ttln' : "ttlnuJet_all22WCsStartPtCheckdim6TopMay20GST_run0",
            'tttt' : "tttt_FourtopsMay3v1_run0",
        }
        self.WCs_per_process = eval( open(os.environ['CMSSW_BASE'] + '/src/CMGTools/TTHAnalysis/data/global_eft/selectedWCs.txt').read())
        self.mods=[]
        self.tmpdirs=[]
        #self.tmpdir='/scratch/'
        self.tmpdir='/scratch365/byates2'
        self.tmpdir=''

        path=os.environ['CMSSW_BASE'] + '/src/CMGTools/TTHAnalysis/data/global_eft/%s/SubProcesses'%self.process_dict[process]


        # reweighting points (first should be reference)
        self.param_cards=[os.environ['CMSSW_BASE'] + '/src/CMGTools/TTHAnalysis/data/global_eft/param_card_%s.dat'%process]
        #self.all_WCs=['cQlMi', 'ctq8', 'ctli', 'cpQM', 'cQq81', 'cQl3i', 'ctlTi', 'cQei', 'ctG', 'ctp', 'cptb', 'cQq13', 'ctZ', 'ctW', 'ctei', 'cpQ3', 'cbW', 'ctt1', 'cQq83', 'cQq11', 'cQQ1', 'cpt', 'ctlSi', 'cQt1', 'cQt8', 'ctq1', 'ctu1', 'cQq1']
        #self.all_WCs=['ctlTi', 'ctq1', 'ctq8', 'clq1', 'cQd8', 'ctu1', 'ctu8', 'cQq83', 'cQQ1', 'cbB', 'cQt1', 'ctb8', 'cbG', 'cld', 'ctt1', 'ctd1', 'cQq81', 'cQlMi', 'cbW', 'cpQ3', 'ctei', 'cQei', 'ctW', 'cpQM', 'ctlSi', 'cQu1', 'ctZ', 'cQl3i', 'ctG', 'clu', 'cQq13', 'cQq11', 'cptb', 'ctli', 'ctd8', 'cQb8', 'ctp', 'cQd1', 'cQu8', 'cQt8', 'cQq18', 'cpt', 'cQq1']
        #self.all_WCs = list(set(self.all_WCs + ["ctp","cpQM","cpQ3","cpt","cptb","ctW","ctZ","cbW","ctG","cQlMi","cQl3i","cQei","ctli","ctei","ctlSi","ctlTi"]))
        self.all_WCs=['ctlTi', 'ctq1', 'ctq8', 'clq1', 'cQd8', 'ctu1', 'ctu8', 'cQq83', 'cQQ1', 'cbB', 'cQt1', 'ctb8', 'cbG', 'cld', 'ctt1', 'ctd1', 'cQq81', 'cQlMi', 'cbW', 'cpQ3', 'ctei', 'cQei', 'ctW', 'cpQM', 'ctlSi', 'cQu1', 'ctZ', 'cQl3i', 'ctG', 'clu', 'cQq13', 'cQq11', 'cptb', 'ctli', 'ctd8', 'cQb8', 'ctp', 'cQd1', 'cQu8', 'cQt8', 'cQq18', 'cpt', 'cQq1']
        self.all_WCs = self.WCs_per_process[process]

        '''
        ttH starting point
        "ctp": -8.000000e-01
        "cpQM": -1.114000e+01
        "cpQ3": 5.790000e+00
        "cpt": 1.288000e+01
        "cptb": 1.653000e+01
        "ctW": 1.000000e+00
        "ctZ": -1.250000e+00
        "cbW": 7.340000e+00
        "ctG": 9.000000e-02
        "cQlMi": 1.000000e+02
        "cQl3i": 1.000000e+02
        "cQei": 1.000000e+02
        "ctli": 1.000000e+02
        "ctei": 1.000000e+02
        "ctlSi": 1.000000e+02
        "ctlTi": 1.000000e+02
        '''

        # generate scan
        param_card_template = open(os.environ['CMSSW_BASE'] + '/src/CMGTools/TTHAnalysis/data/global_eft/param_card_template_%s.dat'%process).read()
        wc_dict = dict( [(wc, '0.000000e+00') for wc in self.all_WCs])
        '''
        wc_dict['ctlTi'] = 100.000000
        wc_dict['ctq1'] = 0.700000
        wc_dict['ctq8'] = 0.930000
        wc_dict['cQlMi'] = 100.000000
        wc_dict['cQq81'] = 0.990000
        wc_dict['cQq83'] = 1.540000
        wc_dict['cbW'] = 7.340000
        wc_dict['cpQ3'] = 5.790000
        wc_dict['ctei'] = 100.000000
        wc_dict['ctlSi'] = 100.000000
        wc_dict['ctW'] = 1.000000
        wc_dict['cpQM'] = -11.140000
        wc_dict['cQei'] = 100.000000
        wc_dict['ctZ'] = -1.250000
        wc_dict['cQl3i'] = 100.000000
        wc_dict['ctG'] = 0.090000
        wc_dict['cQq13'] = 0.680000
        wc_dict['cQq11'] = -0.720000
        wc_dict['cQQ1'] = 0.000000
        wc_dict['cQt1'] = 0.000000
        wc_dict['cQt8'] = 0.000000
        wc_dict['ctt1'] = 0.000000
        wc_dict['ctu1'] = 1.630000
        wc_dict['cptb'] = 16.530000
        wc_dict['ctli'] = 100.000000
        wc_dict['ctp'] = -0.800000
        wc_dict['cpt'] = 12.880000
        wc_dict = dict( [(wc, '0.000000e+00') for wc in wc_dict.keys()])
        '''
        paramcard=param_card_template.format(**wc_dict)
        fo=open('param_card_%s_sm.dat'%(process),'w')
        fo.write(paramcard); fo.close()
        self.param_cards.append('param_card_%s_sm.dat'%(process))


        self.N = len(self.WCs_per_process[process])
        self.names = ['SM'] + self.WCs_per_process[process]
        self.pts = [[1] + list(np.zeros(len(self.WCs_per_process[process])))]
        wc_val = dict( [(wc, "0") for wc in self.all_WCs])
        #self.pts.append([1, float("-8.000000e-01"), float("-1.114000e+01"), float("5.790000e+00"), float("1.288000e+01"), float("1.653000e+01"), float("1.000000e+00"), float("-1.250000e+00"), float("7.340000e+00"), float("9.000000e-02"), float("1.000000e+02"), float("1.000000e+02"), float("1.000000e+02"), float("1.000000e+02"), float("1.000000e+02"), float("1.000000e+02"), float("1.000000e+02")])
        #self.pts_str = ['_'.join([x+'_'+str(int(y)) for x,y in zip(self.WCs_per_process[process], self.pts[-1][1:])])]
        self.pts_str = ['_'.join([x+'_'+y for x,y in wc_val.items()])]
        #self.pts_dict = {'sm': [1] + '_'.join([x+'_'+str(int(y)) for x,y in zip(self.WCs_per_process[process], self.pts[-1][1:])])}
        #self.pts_dict = {'sm': [1] + '_'.join([x+'_'+str(int(y)) for x,y in zip(self.WCs_per_process[process], list(np.zeros(len(self.WCs_per_process[process]))))])}
        n_temp = 0
        n_files = 0
        for i,op in enumerate(self.WCs_per_process[process]):
            if n_files > 0 and n_files % 100 == 0:
                n_temp += 1
                n_files = 0
            if not os.path.isdir(self.tmpdir + 'tmp' + str(n_temp)):
                os.mkdir(self.tmpdir + 'tmp' + str(n_temp))
            wc_dict = dict( [(wc, '0.000000e+00') for wc in self.all_WCs])
            wc_val = dict( [(wc, "0") for wc in self.all_WCs])
            wc_dict[op]="1.000000e+00"
            wc_val[op]="1"
            self.pts.append([1] + [float(x) for wc,x in wc_dict.items() if wc in self.WCs_per_process[process]])
            #self.pts_str.append('_'.join([x+'_'+str(int(y)) for x,y in zip(self.WCs_per_process[process], self.pts[-1][1:])]))
            self.pts_str.append('_'.join([x+'_'+y for x,y in wc_val.items()]))
            paramcard=param_card_template.format(**wc_dict)
            fo=open(self.tmpdir + 'tmp' + str(n_temp) + '/param_card_%s_%s_1.dat'%(process, op),'w')
            fo.write(paramcard); fo.close()
            self.param_cards.append(self.tmpdir + 'tmp' + str(n_temp) + '/param_card_%s_%s_1.dat'%(process, op))
            n_files += 1

            wc_dict = dict( [(wc, '0.000000e+00') for wc in self.all_WCs])
            wc_val = dict( [(wc, "0") for wc in self.all_WCs])
            wc_dict[op]="2.000000e+00"
            wc_val[op]="2"
            self.pts.append([1] + [float(x) for wc,x in wc_dict.items() if wc in self.WCs_per_process[process]])
            #self.pts_str.append('_'.join([x+'_'+str(int(y)) for x,y in zip(self.WCs_per_process[process], self.pts[-1][1:])]))
            self.pts_str.append('_'.join([x+'_'+y for x,y in wc_val.items()]))
            paramcard=param_card_template.format(**wc_dict)
            fo=open(self.tmpdir + 'tmp' + str(n_temp) + '/param_card_%s_%s_2.dat'%(process, op),'w')
            fo.write(paramcard); fo.close()
            self.param_cards.append(self.tmpdir + 'tmp' + str(n_temp) + '/param_card_%s_%s_2.dat'%(process, op))
            n_files += 1

            for i2,op2 in enumerate(self.WCs_per_process[process]):
                if i <= i2: continue
                wc_dict = dict( [(wc, '0.000000e+00') for wc in self.all_WCs])
                wc_val = dict( [(wc, "0") for wc in self.all_WCs])
                wc_dict[op]="1.000000e+00"
                wc_dict[op2]="1.000000e+00"
                wc_val[op]="1"
                wc_val[op2]="1"
                if i != i2: self.pts.append([1] + [float(x) for wc,x in wc_dict.items() if wc in self.WCs_per_process[process]])
                #if i != i2: self.pts_str.append('_'.join([x+'_'+str(int(y)) for x,y in zip(self.WCs_per_process[process], self.pts[-1][1:])]))
                #if i != i2: self.pts_str.append('_'.join([x+'_'+y for x,y in wc_val.items()]))
                self.pts_str.append('_'.join([x+'_'+y for x,y in wc_val.items()]))
                paramcard=param_card_template.format(**wc_dict)
                fo=open(self.tmpdir + 'tmp' + str(n_temp) + '/param_card_%s_%s_%s_1.dat'%(process, op, op2),'w')
                fo.write(paramcard); fo.close()
                self.param_cards.append( self.tmpdir + 'tmp' + str(n_temp) + '/param_card_%s_%s_%s_1.dat'%(process, op, op2) )
                n_files += 1
            '''
            print process
            print self.WCs_per_process[process]
            print len(self.WCs_per_process[process])
            print [wc for wc,x in wc_dict.items() if wc in self.WCs_per_process[process]]
            print len([wc for wc,x in wc_dict.items() if wc in self.WCs_per_process[process]]), "total WCs"
            print len(self.pts[-1]), "points"
            '''


        for card in self.param_cards:
            dirpath = tempfile.mkdtemp(dir=self.tmpdir)
            self.tmpdirs.append(dirpath)
            #print dirpath
            shutil.copyfile( path + '/allmatrix2py.so', self.tmpdirs[-1] + '/allmatrix2py.so')
            sys.path[-1] =self.tmpdirs[-1]
            self.mods.append(imp.load_module('allmatrix2py',*imp.find_module('allmatrix2py')))
            del sys.modules['allmatrix2py']
            print 'initializing', card
            self.mods[-1].initialise(card)
        print self.mods

        self.pdgOrderSorted = [SortPDGs(x.tolist()) for x in self.mods[-1].get_pdg_order()]
        self.pdgOrder = [x.tolist() for x in self.mods[-1].get_pdg_order()]
        self.all_prefix = [''.join(j).strip().lower() for j in self.mods[-1].get_prefix()]
        self.hel_dict = {}; prefix_set = set(self.all_prefix)
        for prefix in prefix_set:
            if hasattr(self.mods[-1], '%sprocess_nhel' % prefix):
                nhel = getattr(self.mods[-1], '%sprocess_nhel' % prefix).nhel
                self.hel_dict[prefix] = {}
                for i, onehel in enumerate(zip(*nhel)):
                    self.hel_dict[prefix][tuple(onehel)] = i + 1


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.wrappedOutputTree = wrappedOutputTree
        for card in self.param_cards:
            self.wrappedOutputTree.branch('weight_%s'%(card.split('/')[-1].replace('param_card_','').replace('.','_')),'F')
        nCols = 1 + 2*self.N + self.N*(self.N - 1)/2
        self.wrappedOutputTree.branch('EFTPostMortem', 'F', nCols)

    def endJob(self):
        for dr in self.tmpdirs:
            shutil.rmtree(dr, ignore_errors=True)

    def analyze(self, event):

        lheParts = [l for l in Collection(event, 'LHEPart')]
        pdgs = [x.pdgId for x in lheParts]
        hel  = [x.spin  for x in lheParts]



        p = [ ]
        for part in lheParts:
            if part.status < 0:
                energy = math.sqrt(part.incomingpz*part.incomingpz+part.mass*part.mass)
                p.append([energy,0.,0.,part.incomingpz])
            else:
                p.append([part.p4().E(), part.p4().Px(), part.p4().Py(), part.p4().Pz()])

        # madgraph pads processes with low multiplicities
        while len(pdgs) < len(self.pdgOrderSorted[0]):
            pdgs.append(0)

        evt_sorted_pdgs = SortPDGs(pdgs)

        try:
            idx = self.pdgOrderSorted.index(evt_sorted_pdgs)
        except ValueError:
            print self.pdgOrderSorted
            print '>> Event with PDGs %s does not match any known process' % pdgs
            return res

        target_pdgs=self.pdgOrder[idx]
        pdgs_withIndices = [(y,x) for x,y in enumerate(pdgs)]
        mapping=[]


        for p1 in target_pdgs:
            toremove=None
            for p2 in pdgs_withIndices:
                if p2[0]==p1:
                    mapping.append( p2[1])
                    toremove=p2
                    break
            if toremove:
                pdgs_withIndices.remove(toremove)
            else:
                raise RuntimeError("It shouldn't be here")

        final_pdgs = []
        final_parts = []
        final_hels = []
        for in_Indx in mapping:
            final_pdgs.append(pdgs[in_Indx])
            if final_pdgs[-1] == 0: continue
            final_parts.append(p[in_Indx])
            final_hels.append(hel[in_Indx])

        if target_pdgs != final_pdgs:
            raise RuntimeError("Wrong pdgid")

        hel_dict = self.hel_dict[self.all_prefix[idx]]
        t_final_hels = tuple(final_hels)

        if t_final_hels in hel_dict:
            nhel = hel_dict[t_final_hels]
        else:
            print "Available helicities are"
            print hel_dict
            print "tried", t_final_hels
            raise RuntimeError("Helicity configuration not found")

        com_final_parts = []


        pboost = [final_parts[0][i] + final_parts[1][i] for i in xrange(4)]

        for part in final_parts:
            com_final_parts.append(zboost(part, pboost))


        final_parts_i = invert_momenta(com_final_parts)
        scale2=0
        weights=[]
        for mod in self.mods:
            if 0 in final_pdgs: final_pdgs.remove(0)
            weights.append( mod.smatrixhel( final_pdgs, final_parts_i, event.LHE_AlphaS, scale2, nhel) )
            #print 'Adding weight:', weights[-1]/weights[0]

        nCols = 1 + 2*self.N + self.N*(self.N - 1)/2
        nRows = len(self.pts)
        pairs = []
        for row_idx in range(len(self.pts[0])):
            for col_idx in range(len(self.pts[0])):
                if col_idx > row_idx: continue
                pairs.append((row_idx, col_idx))

        '''
        import numpy as np

        print(nRows, nCols, pairs)

        # Parse all pts_str once up front to avoid doing it per-row
        parsed_pts = [
            {k: int(v) for k, v in zip(x.split('_')[::2], x.split('_')[1::2])}
            for x in self.pts_str
        ]

        # Prepare A and b
        A = np.zeros((nRows, nCols), dtype=float)
        b = np.zeros(nRows, dtype=float)

        # Extract names from pairs
        n1_list = [self.names[pair[0]] for pair in pairs]
        n2_list = [self.names[pair[1]] for pair in pairs]

        for row_idx in range(nRows):
            d = parsed_pts[row_idx]
            x1_vals = np.array([1.0 if n1 == 'SM' else d[n1] for n1 in n1_list], dtype=float)
            x2_vals = np.array([1.0 if n2 == 'SM' else d[n2] for n2 in n2_list], dtype=float)
            A[row_idx, :] = x1_vals * x2_vals
            b[row_idx] = weights[row_idx + 1] / weights[0]
        '''

        A = np.zeros((nRows, nCols))
        b = np.zeros(nRows)
        '''
        print nRows, nCols, pairs
        for row_idx in range(nRows):
            for col_idx in range(nCols):
                idx_pair = pairs[col_idx]
                n1 = self.names[idx_pair[0]]
                n2 = self.names[idx_pair[1]]
                #print n1, n2, 'idx_pair=', idx_pair, 'row_idx=', row_idx
                #print self.N, 'WCs'
                #print self.pts[row_idx]
                #print len(self.pts[row_idx])
                #print len(self.pts)
                x1 = 1.0 if n1 == 'SM' else self.pts[row_idx][idx_pair[0]]# Hard set SM value to 1.0
                x2 = 1.0 if n2 == 'SM' else self.pts[row_idx][idx_pair[1]]# Hard set SM value to 1.0
                #print(n1, x1, n2, x2)
                #print [dict(zip(x.split('_')[::2], x.split('_')[1::2]))   for x in self.pts_str]
                #d = [dict(zip(x.split('_')[::2], x.split('_')[1::2]))   for x in self.pts_str][row_idx]
                #x1 = 1.0 if n1 == 'SM' else int(d[n1])# Hard set SM value to 1.0
                #x2 = 1.0 if n2 == 'SM' else int(d[n2])# Hard set SM value to 1.0
                #print(n1, x1, n2, x2)
                #d = [dict(zip(x.split('_')[::2], x.split('_')[1::2]))   for x in self.pts_str][row_idx]
                #x1 = 1.0 if n1 == 'SM' else int(d[n1])# Hard set SM value to 1.0
                #x2 = 1.0 if n2 == 'SM' else int(d[n2])# Hard set SM value to 1.0

                A[row_idx,col_idx] = x1*x2
                #b[row_idx] = weights[row_idx+1]#*weights[1]/weights[0]
                #b[row_idx] = weights[row_idx+1]/weights[1]*weights[1]/weights[0];
                b[row_idx] = weights[row_idx+1]/weights[0];
        '''
        parsed_pts = [
            {k: float(v) for k, v in zip(x.split('_')[::2], x.split('_')[1::2])}
            for x in self.pts_str
        ]

        # Extract names from pairs
        n1_list = [self.names[pair[0]] for pair in pairs]
        n2_list = [self.names[pair[1]] for pair in pairs]

        LHEWeight_originalXWGTUP = 0.0069839

        for row_idx in range(nRows):
            #row = np.asarray(self.pts[row_idx])
            d = parsed_pts[row_idx]
            x1_vals = np.array([1.0 if n1 == 'SM' else d[n1] for n1 in n1_list], dtype=float)
            x2_vals = np.array([1.0 if n2 == 'SM' else d[n2] for n2 in n2_list], dtype=float)
            A[row_idx, :] = x1_vals * x2_vals
            b[row_idx] = weights[row_idx + 1] / weights[0]# * LHEWeight_originalXWGTUP

        #print nRows, nCols, pairs
        '''
        print(self.pts_str)
        print(self.param_cards)
        print(A, b)
        exit()
        '''

        for i, card in enumerate(self.param_cards):
            self.wrappedOutputTree.fillBranch('weight_%s'%(card.split('/')[-1].replace('param_card_','')).replace('.','_'), weights[i]/weights[0])

        x, resid, _, _ = np.linalg.lstsq(A,b,rcond=1e-20)
        #resid = b - np.dot(A, x)
        #print(resid)
        # Testing with ROOT, results were worse
        #buf = array.array('d', A.flatten())
        #B = TMatrixD(nRows, nCols, buf)
        #svd = TDecompSVD(B)
        #bvec = TVectorD(len(b), array.array('d', b))
        #ok = svd.Solve(bvec)
        #x = np.array([bvec[i] for i in range(bvec.GetNrows())])
        #x, resid, _, _ = np.linalg.lstsq(A,b,rcond=None)
        '''
        print('A', A)
        print('b', b)
        print('x', x)
        print('A.x', A.dot(x))
        print(A.dtype, b.dtype, x.dtype)
        exit()
        '''
        self.wrappedOutputTree.fillBranch('EFTPostMortem', x.T)


        return True

eft_TTll = lambda : globalEFTreWeighting('ttll')
eft_TTH  = lambda : globalEFTreWeighting('ttH')
eft_TTTT = lambda : globalEFTreWeighting('tttt')
eft_THQ  = lambda : globalEFTreWeighting('tHq')
eft_TllQ = lambda : globalEFTreWeighting('tllq')
eft_TTln = lambda : globalEFTreWeighting('ttln')
