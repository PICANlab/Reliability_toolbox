#!/bin/csh
# does all cbt single subject regressions with area under the curve modelling for the SRET task
# assumes we have regressors made 
# This script assumes we have the AFNI in our path
# This code has been described in in Compere et al (2020).
# note for 3dDeconvolve censor files 0=censored. 1=use.


set AFNIPATH = '/etc/afni/abin'

set fdir = /data/Siegle/reliability/ssanal/sret/ssregs_canonical_amplitude
mkdir -p $fdir
cd $fdir

set tasks = (cbt cbtbirc)

foreach task ($tasks)
foreach sub (`ls /data/Siegle/${task}/dfs/ | grep sretnegtrials.reg | cut -c 1-4`)
 if (-e ${fdir}/sret${sub}.ssreg+orig.BRIK) then
   echo using existing ssreg sret for $sub
 else

   if (-e /data/Siegle/cbt/grptrans/sret/cbt${sub}/func_smooth/grptranssret+orig.HEAD) then
     set inp = /data/Siegle/cbt/grptrans/sret/cbt${sub}/func_smooth/grptranssret+orig.HEAD
     set pre = PRES_CBT
     set this_task = cbt
   else
     if (-e /data/Siegle/cbtbirc/grptrans/sret/cbt${sub}/func_smooth/grptranssret+orig.HEAD) then 
         set inp = /data/Siegle/cbtbirc/grptrans/sret/cbt${sub}/func_smooth/grptranssret+orig.HEAD
         set pre = BIRC_CBT
         set this_task = cbtbirc
     else
       set inp = NONEFOUND
       set pre = NOWHERE
       echo "No Group Trans Data Found for sret $sub anywhere (cbt, cbtbirc)"
     endif
   endif
        
   if (-e $inp) then
     echo ==========================================
     echo WORKING ON SSREG SRET $sub $pre NEGATIVE TRIALS
     cd $fdir
     ${AFNIPATH}/3dDeconvolve \
      -input $inp \
      -mask /data/Siegle/cbt/struct/refbrain/minicolin+orig \
      -num_stimts 1 \
      -nfirst 1 \
      -stim_times 1 ../../../stim_times/${sub}_stim_times.1D 'BLOCK5(1,1)' -stim_label 1 negtrials \
      -censor "/data/Siegle/${this_task}/dfs/${sub}sretnegtrials.cenafni" \
      -polort 0 \
      -fout -rout -tout \
      -full_first \
      -bucket ${fdir}/sret${sub}${pre}-negtrials.ssreg \
      -GOFORIT 2
    endif
  endif
end
end
