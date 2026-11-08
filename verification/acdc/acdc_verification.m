% KLUSTR v1.0: ACDC STEADY STATES FOR THE VERIFICATION (Sect. 5 of the paper) by R. Halonen 2026
% Sulfuric acid (A) and dimethylamine (D), cluster sets 4x4 (input_AD.inp) and band8 (input_AD_band8.inp),
% free energies olenius2013 (H and S) and myllys19_anharmonic (G at each temperature).
% Acid monomer (1A + 1A1D + ... + 1A4D) and free base fixed, sink CS (d/d_1A)^-1.7.
%
% Output acdc_<free energies>_<set>.csv, one row per condition: T, CS, CA, CB, converged (1 yes, 0 no, -1 negative
% concentrations), J (cm-3 s-1), time (s), stable (1: the largest cluster collides more often than it evaporates)
% and the concentrations of all clusters (cm-3).

Ts=208:10:328;             % K
CAs=10.^(3:0.5:11);        % acid monomer, cm-3
CBs=10.^(2:0.5:12);        % free base, cm-3
CSs=[0 10.^(-5:0.5:0)];    % sink of 1A, s-1
tmax=120;                  % s, wall-time limit of one steady state

cd(fileparts(mfilename('fullpath')))
fclose('all');
% steady state: windows of 20 output times, converged when no concentration changes by more than 1e-6
opts={'Sources_in','sources.txt','difftol',1e-6,'Options','''RelTol'',1e-8,''AbsTol'',1e-4,''Events'',@time_limit'};
tw=@(L) L*(1:20)/20;

for fe=["olenius2013" "myllys19_anharmonic"]
    for cset=["4x4" "band8"]
        if cset=="4x4", inp='input_AD.inp'; big="4A4D"; else, inp='input_AD_band8.inp'; big="8A8D"; end
        fid=-1; times=[]; nbad=0; nstop=0; t1=tic;
        for T=Ts
            tT=tic;
            if fe=="olenius2013", energies='SA-DMA_olenius2013_HS.txt'; else, energies=sprintf('%s/G%dK.txt',fe,T); end
            [status,msg]=system(sprintf(['perl ACDC/acdc_2024_02_12.pl --i %s --e %s --temperature %d --cs exp_loss ' ...
                '--exp_loss_coefficient 1 --exp_loss_exponent -1.7 --silent'],inp,energies,T));
            if status, error('ACDC generator failed at %d K:\n%s',T,msg), end
            clear driver_acdc equations_acdc dofluxes formationrate get_coll get_evap get_cs get_wl
            rehash

            % the largest cluster (the last one in ACDC's order): collision coefficients with every cluster
            % (cm3 s-1) and total evaporation rate (s-1), every channel once as in get_evap.m
            get_coll
            Kb=1e6*K(end,:); Kb(isnan(Kb))=0;
            L=readlines('get_evap.m'); Eb=0;
            for l=L(contains(L,"% "+big+" ->"))'
                Eb=Eb+str2num(char(extractBetween(l,'= ',';'))); %#ok<ST2NM>
            end

            for CS=CSs
                for CB=CBs
                    for CA=CAs
                        writelines(sprintf('constant 1A %e -1A1D-1A2D-1A3D-1A4D\nconstant 1D %e',CA,CB),'sources.txt');
                        t0=tic;
                        time_limit('set',tmax);
                        % from no clusters, then longer and longer windows until converged
                        [ok,cl,C]=steady(tw(1e-8),[],opts,CS,{});
                        J=NaN;
                        if ok~=-1, [ok,cl,C,J]=steady(tw(1e-1),C,opts,CS,{}); end
                        if ok~=-1
                            for e=[4 4 4 3*ones(1,20) 5 6 7 8]
                                [ok,cl,C,J]=steady(tw(10^e),C,opts,CS,{'repeat'});
                                if ok==1 || ok==-1 || time_limit('expired'), break, end
                            end
                        end
                        stopped=time_limit('expired');
                        if ok==-1, J=NaN; elseif stopped, ok=0; end
                        time=toc(t0);
                        if ~isempty(cl), clust=cl; end
                        if fid<0
                            fid=fopen(sprintf('acdc_%s_%s.csv',fe,cset),'w');
                            fprintf(fid,'T,CS,CA,CB,converged,J,time,stable%s\n',sprintf(',%s',clust{:}));
                        end
                        n=numel(clust);
                        if numel(C)<n, C=nan(1,n); end
                        stable=double(Eb<=Kb(1:n)*reshape(C(1:n),[],1));
                        fprintf(fid,'%d,%g,%e,%e,%d,%e,%.4f,%d%s\n',T,CS,CA,CB,ok,J(end),time,stable,sprintf(',%e',C(1:n)));
                        times(end+1)=time; nbad=nbad+(ok~=1); nstop=nstop+stopped; %#ok<SAGROW>
                        if mod(numel(times),200)==0
                            fprintf('%s %s %d K: %d conditions, %d not converged, %.0f s\n',fe,cset,T,numel(times),nbad,toc(tT));
                        end
                    end
                end
            end
        end
        fclose(fid);
        fprintf('%s %s: %d conditions, %d not converged, %d stopped at %d s, time per steady state %.3f s (mean), %.3f s (median), %.1f h\n', ...
            fe,cset,numel(times),nbad,nstop,tmax,mean(times),median(times),toc(t1)/3600);
    end
end
delete('driver_acdc.m','equations_acdc.m','dofluxes.m','formationrate.m','get_*.m','sources.txt')


function [ok,clust,C,J]=steady(tw,C0,opts,CS,extra)
% One call of ACDC's driver. With negative concentrations the driver returns early without J: ok = -1.
    clust={}; J=NaN;
    try
        if isempty(C0)
            [~,~,ok,clust,C,~,~,J]=driver_acdc(tw,opts{:},'CSfactor',CS,extra{:});
        else
            [~,~,ok,clust,C,~,~,J]=driver_acdc(tw,C0,opts{:},'CSfactor',CS,extra{:});
        end
    catch err
        if ~strcmp(err.identifier,'MATLAB:unassignedOutputs'), rethrow(err), end
        ok=-1; C=nan(size(C0));
    end
end
