# Generates repro_old.do and repro_new.do (same config blocks, mirroring test/suite/mkrun.py's
# per-config block verbatim; the new file adds two fullprecisionslow legs as discriminators).
import os, sys
_REPO = os.environ.get("WOOLDID_REPO", "")
_W    = os.environ.get("WOOLDID_SCRATCH", "")
if not _REPO: raise SystemExit("set WOOLDID_REPO=<path>  (see repro/README.md)")
if not _W: raise SystemExit("set WOOLDID_SCRATCH=<path>  (see repro/README.md)")
W=_W
R=_REPO
DTA_CONT=f'{R}/test/data/dta_cont.dta'
DTA_GAP=f'{R}/test/data/dta_path_relyrgap.dta'
DTA_RCS=f'{R}/test/suite/data/survey_rcs.dta'

def block(tag, cfg, dataset, weight, prep, options, prints, sentinels, dumpdir, extra_opt=''):
    L=[f'* ======== {tag}: config {cfg} ========','fresh']
    L.append(f'use "{dataset}", clear' if dataset else 'clear')
    L.extend(prep)
    call='cap noi wooldid y i t ttre' if cfg!='ext_contract10_rcs_factor_fe' else 'cap noi wooldid y_rcs state_group t ttre'
    if weight: call+=f' [{weight}]'
    call+=', '+options+((' '+extra_opt) if extra_opt else '')
    L.append(call)
    L.append('local rc = _rc')
    L.append(f'dis "{tag}_rc=`rc\'"')
    L.append(f'dis "{tag}_cmdline=[`e(cmdline)\']"')
    L.append(f'dis "{tag}_properties=[`e(properties)\']"')
    dumpname = cfg if not extra_opt else cfg+'_fps'
    L.append(f'dumpres {dumpname} `rc\' "{dumpdir}"')
    for p in prints: L.append(p.replace('TAG',tag))
    for s in sentinels: L.append(f'dis "{s}"')
    L.append('')
    return '\n'.join(L)

def prints_cont_sg():
    return ['pcell TAG_main_ame_est wdidmainresults 3 1','pcell TAG_main_ame_se wdidmainresults 3 2',
            'pcell TAG_main_ame_t wdidmainresults 3 3','pcell TAG_main_ame_p wdidmainresults 3 4',
            'pcell TAG_ESame_pre2_se wdidmainresults 9 2','pcell TAG_ESame_post1_se wdidmainresults 12 2',
            'pcell TAG_main_att_se wdidmainresults 1 2',
            'pcell TAG_V_3_3 V 3 3','pcell TAG_V_9_9 V 9 9','pcell TAG_V_1_1 V 1 1',
            'pcell TAG_pwsgs_ame_est pwsgs 3 1','pcell TAG_pwsgs_ame_se pwsgs 3 2','pcell TAG_pwsgs_ame_t pwsgs 3 3','pcell TAG_pwsgs_ame_p pwsgs 3 4',
            'pcell TAG_pwsgs_att_se pwsgs 1 2',
            'pscal TAG_joint_Choriz joint_Choriz','pscal TAG_joint_horiz joint_horiz',
            'pscal TAG_joint_horiz_s1 joint_horiz_s1','pscal TAG_joint_horiz_s2 joint_horiz_s2',
            'pscal TAG_joint_pre_att_allsgssame joint_pre_att_allsgssame','pscal TAG_joint_ame_allsgssame joint_ame_allsgssame',
            'pscal TAG_wdid_fcol wdid_fcol']
def prints_cont_lat():
    return ['pscal TAG_rss rss','pscal TAG_rmse rmse','pscal TAG_r2 r2','pscal TAG_r2adj r2adj','pscal TAG_r2_within r2_within',
            'pscal TAG_r2_withinadj r2_withinadj','pscal TAG_ll ll','pscal TAG_ll_0 ll_0','pscal TAG_df_r df_r','pscal TAG_N N',
            'pname TAG_main_ame main_ame','pname TAG_main_att main_att','pname TAG_ES_ame_pre2 ES_ame_pre2',
            'pcell TAG_main_att_se wdidmainresults 1 2','pcell TAG_main_ame_se wdidmainresults 3 2','pscal TAG_wdid_fcol wdid_fcol']
def prints_d3():
    return ['pcell TAG_main_att_est wdidmainresults 1 1','pcell TAG_main_att_se wdidmainresults 1 2',
            'pcell TAG_main_ame_est wdidmainresults 3 1','pcell TAG_main_ame_se wdidmainresults 3 2',
            'pcell TAG_ESame_pre2_est wdidmainresults 9 1','pcell TAG_ESame_post1_est wdidmainresults 12 1',
            'pname TAG_b_main_ame main_ame','pname TAG_b_main_att main_att','pcell TAG_V_3_3 V 3 3','pcell TAG_V_1_1 V 1 1']
def prints_gap():
    return ['pname TAG_b_main_SG1_att main_SG1_att','pname TAG_b_main_SG1_ame main_SG1_ame','pname TAG_b_pretreatment_SG0_att pretreatment_SG0_att',
            'pname TAG_b_ES_ame_SG2_post1 ES_ame_SG2_post1','pname TAG_b_ES_att_SG0_post1 ES_att_SG0_post1',
            'pshape TAG_sg0_shape wdidsgresults_0','pshape TAG_sg1_shape wdidsgresults_1',
            'pcell TAG_sg0_r2_est wdidsgresults_0 2 1','pcell TAG_sg0_r2_relyr wdidsgresults_0 2 7','pcell TAG_sg0_r3_est wdidsgresults_0 3 1',
            'pcell TAG_sg0_r8_est wdidsgresults_0 8 1',
            'pcell TAG_pwsgs_att_est pwsgs 1 1','pcell TAG_pwsgs_att_se pwsgs 1 2','pcell TAG_pwsgs_att_t pwsgs 1 3','pcell TAG_pwsgs_att_p pwsgs 1 4',
            'pcell TAG_pwsgs_ame_est pwsgs 7 1','pcell TAG_pwsgs_ame_se pwsgs 7 2',
            'pcell TAG_main_ame_est wdidmainresults 3 1','pcell TAG_main_att_se wdidmainresults 1 2',
            'pcell TAG_V_1_1 V 1 1','pcell TAG_V_3_3 V 3 3',
            'pscal TAG_joint_att_allsgssame joint_att_allsgssame','pscal TAG_wdid_sgdegen wdid_sgdegen','pscal TAG_wdid_sgdegenr wdid_sgdegenr',
            'pscal TAG_wdid_fcguard wdid_fcguard','pscal TAG_wdid_fcol wdid_fcol']
def prints_rcs():
    return ['pname TAG_b_main_att main_att','pname TAG_b_pretreatment_SG1_att pretreatment_SG1_att','pname TAG_b_ES_att_SG1_pre2 ES_att_SG1_pre2',
            'pname TAG_b_ES_att_SG1_post0 ES_att_SG1_post0','pname TAG_b_ES_att_SG2_pre2 ES_att_SG2_pre2',
            'pshape TAG_sg1_shape wdidsgresults_1','pshape TAG_sg2_shape wdidsgresults_2',
            'pcell TAG_sg1_r2_est wdidsgresults_1 2 1','pcell TAG_sg1_r2_relyr wdidsgresults_1 2 7','pcell TAG_sg1_r3_est wdidsgresults_1 3 1',
            'pcell TAG_main_att_se wdidsgresults_1 1 2','pcell TAG_pooled_main_att_se wdidmainresults 1 2',
            'pscal TAG_wdid_fcguard wdid_fcguard','pscal TAG_wdid_fcol wdid_fcol','pscal TAG_N_clust N_clust','pscal TAG_N N']

def header(builddir):
    return f'''* Reproduction, contsg group (stage B proofs). Generated by gen_repro.py - edit that file, not this one.
* Mirrors test/suite/mkrun.py's per-config block (seed 98765 / sortseed 98765 / ereturn clear / drops / use / preps.json lines).
clear all
set more off
set maxvar 20000
set graphics off
adopath ++ "{builddir}"
include "{W}/repro_helpers.doh"
dis "BUILD_DIR={builddir}"
which wooldid
'''

def body(dumpdir, extras):
    B=[]
    B.append(block('AMESE_R1','cont_sg',DTA_CONT,'',['qui gen sg = mod(i,2)+1'],
        'contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests safety(off)',
        prints_cont_sg(),['REPRO_AMESE_R1_DONE'] + ([] if extras else ['REPRO_QRSGJT_DONE']),dumpdir))
    B.append(block('R4W','cont_lat',DTA_CONT,'aw=w',[],
        'contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2)',
        prints_cont_lat(),['REPRO_R4W_DONE'],dumpdir))
    B.append(block('SIA4','trig_d3void',None,'',
        ['qui set obs 12','qui gen i = _n','qui gen ttre = cond(_n <= 4, 5, cond(_n <= 8, 7, .))','qui expand 4','qui bysort i: gen unit = _n',
         'qui gen dose = 0.5 + 0.1*i','qui expand 10','qui bysort i unit: gen t = _n',
         'qui gen y = 0.05*i + 0.04*t + 0.3*dose*(ttre < . & t >= ttre) + 0.3*sin(2.2*i + 1.1*t + 0.7*unit)','sort i unit t','ereturn clear'],
        'att robust contreat(dose) contreatcontrols(basic) esprelength(2) espostlength(2)',
        prints_d3(),['REPRO_SIA4_DONE'],dumpdir))
    B.append(block('SGLAB','cont_sg_gap',DTA_GAP,'',['gen dose2 = 0.5 + mod(i,5) + 0.2*mod(i,3)','gen sgr = mod(i,3)'],
        'contreat(dose2) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sgr)',
        prints_gap(),['REPRO_SGLAB_DONE','REPRO_BS8_DONE','REPRO_A17SGV_DONE','REPRO_SIA4GP_DONE'],dumpdir))
    B.append(block('SI7','ext_contract10_rcs_factor_fe',DTA_RCS,'pw=survey_weight',['replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .'],
        'att subgroups(group) fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) safety(off)',
        prints_rcs(),['REPRO_SI7_DONE'] if not extras else [],dumpdir))
    if extras:
        B.append(block('SI7FPS','ext_contract10_rcs_factor_fe',DTA_RCS,'pw=survey_weight',['replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .'],
            'att subgroups(group) fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) safety(off)',
            prints_rcs(),['REPRO_SI7_DONE'],dumpdir+'/extra',extra_opt='fullprecisionslow'))
        B.append(block('QRFPS','cont_sg',DTA_CONT,'',['qui gen sg = mod(i,2)+1'],
            'contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests safety(off)',
            ['pscal TAG_joint_horiz_s1 joint_horiz_s1','pscal TAG_joint_horiz_s2 joint_horiz_s2',
             'pscal TAG_joint_pre_att_allsgssame joint_pre_att_allsgssame','pscal TAG_joint_horiz joint_horiz','pscal TAG_wdid_fcol wdid_fcol',
             'pcell TAG_main_att_se wdidmainresults 1 2'],
            ['REPRO_QRSGJT_DONE'],dumpdir+'/extra',extra_opt='fullprecisionslow'))
    return '\n'.join(B)

open(f'{W}/repro_old.do','w').write(header(f'{W}/old')+body(f'{W}/old_dumps',False)+'dis "REPRO_ALL_OLD_DONE"\nexit, clear\n')
open(f'{W}/repro_new.do','w').write(header(f'{W}/new')+body(f'{W}/new_dumps',True)+'dis "REPRO_ALL_NEW_DONE"\nexit, clear\n')
print('written')
