% KLUSTR v1.0 verification (Sect. 5 of the paper): ACDC steady states for sulfuric acid (A) and dimethylamine (D),
% cluster sets 4x4 (input_AD.inp) and band8 (input_AD_band8.inp), free energies olenius2013 and myllys19_anharmonic.
% Acid monomer (1A + 1A1D + ... + 1A4D) and free base fixed, sink CS (d/d_1A)^-1.7. Run in this folder: acdc_verification
% Writes acdc_<free energies>_<set>.csv: T, CS, CA, CB, converged, J, time (s, of the steady state), stable (1: the
% largest cluster of the set collides more often than it evaporates in the steady state) and all cluster concentrations.

energy_sets={'olenius2013','myllys19_anharmonic'};
sets={'4x4','band8'};
Ts=208:10:328;             % K
CAs=10.^(3:0.5:11);        % acid monomer, cm-3
CBs=10.^(2:0.5:12);        % free base, cm-3
CSs=[0 10.^(-5:0.5:0)];    % sink of 1A, s-1
tmax=120;                  % s, wall-time limit of one steady state

fclose('all');
home=fileparts(mfilename('fullpath'));   % the folder of this script, wherever MATLAB is
acdc=fullfile(home,'ACDC','acdc_2024_02_12.pl');   % ACDC version 2024_02_12 (GPL-3.0, ACDC/LICENSE)
cd(home)
addpath(home)              % time_limit.m, also from the work folders

% steady-state settings: windows of 20 output times, the driver's convergence test at 1e-6
opts={'Sources_in','sources.txt','difftol',1e-6,'Options','''RelTol'',1e-8,''AbsTol'',1e-4,''Events'',@time_limit'};
tw=@(L) L*(1:20)/20;

for fe=energy_sets
    fe=fe{1};
    for set=sets
        set=set{1};
        inp=fullfile(home,'input_AD.inp');
        if strcmp(set,'band8'), inp=fullfile(home,'input_AD_band8.inp'); end
        out=fullfile(home,['acdc_' fe '_' set '.csv']);

        % ACDC writes its files into the current folder, so each run works in its own folder
        work=fullfile(home,['work_' fe '_' set]);
        if ~exist(work,'dir'), mkdir(work); end
        cd(work)

        fid=[]; clust={}; times=[]; nbad_all=0; nstop=0; tset=tic;
        for T=Ts
            tT=tic;
            if strcmp(fe,'olenius2013')
                energies=fullfile(home,'SA-DMA_olenius2013_HS.txt');
            else
                energies=fullfile(home,fe,sprintf('G%dK.txt',T));
            end
            % equations for this temperature; the sink is generated for CS = 1 s-1 and scaled with CSfactor
            [status,msg]=system(sprintf(['perl "' acdc '" --i "%s" --e "%s" --temperature %d --cs exp_loss --exp_loss_coefficient 1 ' ...
                '--exp_loss_exponent -1.7 --silent'],inp,energies,T));
            if status~=0, cd(home), error('ACDC generator failed at %d K:\n%s',T,msg); end
            clear driver_acdc equations_acdc dofluxes formationrate get_coll get_evap get_cs get_wl
            rehash
            [K,E]=coefficients(); Eb=[];
            k=0; nbad=0; n=numel(CSs)*numel(CBs)*numel(CAs);
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
                        if isempty(fid) && ~isempty(clust)
                            fid=fopen(out,'w');
                            fprintf(fid,'T,CS,CA,CB,converged,J,time,stable%s\n',sprintf(',%s',clust{:}));
                        end
                        if numel(C)<numel(clust), C=nan(1,numel(clust)); end
                        if isempty(Eb) && ~isempty(clust), [Eb,Kb]=largest(K,E,clust); end
                        % the largest cluster: evaporation against collisions with every cluster (s-1)
                        stable=double(all(isfinite(C(1:numel(clust)))) && Eb<=Kb(1:numel(clust))*reshape(C(1:numel(clust)),[],1));
                        fprintf(fid,'%d,%g,%e,%e,%d,%e,%.4f,%d%s\n',T,CS,CA,CB,ok,J(end),time,stable,sprintf(',%e',C(1:numel(clust))));
                        k=k+1; nbad=nbad+(ok~=1); nstop=nstop+stopped; times(end+1)=time; %#ok<SAGROW>
                        if mod(k,200)==0
                            fprintf('%s %s %d K: %d/%d conditions, %d not converged, %.0f s\n',fe,set,T,k,n,nbad,toc(tT));
                        end
                    end
                end
            end
            fprintf('%s %s %d K done: %d conditions, %d not converged, %.0f s\n',fe,set,T,n,nbad,toc(tT));
            nbad_all=nbad_all+nbad;
        end
        fclose(fid);
        fprintf(['%s %s: %d conditions, %d not converged, %d stopped at %d s; time per steady state: mean %.3f s, ' ...
            'median %.3f s; total %.2f h\n'],fe,set,numel(times),nbad_all,nstop,tmax,mean(times),median(times),toc(tset)/3600);
        cd(home)
    end
end


function [K,E]=coefficients()
% ACDC's collision (K, m3 s-1) and evaporation (E, s-1) coefficients, from the files of the generator in this folder
    get_coll
    get_evap
end


function [Eb,Kb]=largest(K,E,clust)
% The largest cluster of the set: its total evaporation rate Eb (s-1), each channel once as in get_evap.m, and its
% collision coefficients with every cluster Kb (cm3 s-1; 0 for collisions that are not in the set)
    v=zeros(numel(clust),2);
    for i=1:numel(clust)
        t=regexp(clust{i},'(\d+)A','tokens','once'); if ~isempty(t), v(i,1)=str2double(t{1}); end
        t=regexp(clust{i},'(\d+)D','tokens','once'); if ~isempty(t), v(i,2)=str2double(t{1}); end
    end
    [~,b]=max(sum(v,2));
    Eb=0;
    for i=1:numel(clust)
        for j=i:numel(clust)
            if all(v(i,:)+v(j,:)==v(b,:)) && isfinite(E(i,j)), Eb=Eb+E(i,j); end
        end
    end
    Kb=1e6*K(b,:); Kb(~isfinite(Kb))=0;
end


function [ok,clust,C,J]=steady(tw,C0,opts,CS,extra)
% One call of ACDC's driver. When the driver finds negative concentrations it returns early without
% the formation rate; the condition is then marked with ok = -1.
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
