% Matlab class to implement a simple retirement model with consumption and savings,
% income shocks and credit constraint
% Written by Fedor Iskhakov, Australian National University, 2016
% See Iskhakov, Jorgensen, Rust and Schjerning 
% "The Endogenous Grid Method for Discrete-Continuous Dynamic Choice Models 
%  with (or without) Taste Shocks" (QE, 2017)


classdef my_model_migration < handle
% This class defines a Deaton consumption-savings model

properties (Access=public)
	%Default parameters of the model
	label		= 'Consumption model with migration decisions'; %name of this model
	Tbar		= 25			; %number of periods (fist period is t=1) 
	ngridm	= 500		; %number of grid points over assets
	mmax		= 50		; %maximum level of assets
	expn		=	5 		; %number of quadrature points used in calculation of expectations
	nsims		= 10		; %number of simulations
	init    =[10 30] ; %interval of the initial wealth
	r 			= 0.05	; %interest rate
	df 			= 0.95	; %discount factor
	sigma   = 0.25	; %sigma parameter in income shocks
	%dus     =	0.005	; %disutility of staying
	theta		= 1.95 	; %CRRA coefficient (log utility if ==1)
	inc0		= 0.75	; %income equation: constant
	inc1		= 0.04  ; %income equation: age coef
	inc2		= 0.0002;	%income equation: age^2 coef
    um0        = 0.005;
    um1        = - 0.0005;
    %um2        = - 0.0005
	cfloor  =	0.001	; %consumption floor (safety net in retirement)
	lambda  = 0.2 ; %scale of the EV taste shocks 
    Tr = 22 ;    %retirement age
end %properties

properties (SetAccess=private, GetAccess=public)
	%Entities computed inside of the model
	policy=polyline; %optimal consumption policy
	value=polyline;
	sims=struct();
end %properties

% working = 1 for working, 2 for retiree

methods (Access=public)
	%Definition of the model
	%Vectorize!!!
    function um=mutility(me,it) %income in period it with given normal shock
		%assume it=1 is age=20
		age=it+19;
        um=me.um0 + me.um1*age %+ me.um2*age*age;
	end %income
    function u=util(me,consumption,choice,it) %utility
		if me.theta==1
			u=log(consumption);
		else
			u=(consumption.^(1-me.theta)-1)/(1-me.theta);
		end
		u=u + mutility(me,it)*(choice==2);
        %u=u - me.dus*(choice==1);
	end %util
	function mu=mutil(me,consumption) %marginal utility
		if me.theta==1
			mu=1./consumption;
		else
			mu= consumption.^(-me.theta);
		end
	end %mutil
	function cons=imutil(me,mutil) %inverse marginal utility
		if me.theta==1
			cons=1./mutil;
		else
			cons= mutil.^(-1/me.theta);
		end
	end %imutil
	function w=income(me,it,shock) %income in period it with given normal shock
		%assume it=1 is age=20
		age=it+19;
        if it < me.Tr 
            w=exp(me.inc0 + me.inc1*age - me.inc2*age.*age + shock);
        else 
            w= 1
        end
	end %income
    function w1=budget(me,it,savings,shocks,choice) 
		%wealth in period t+1, where it=t
		%inputs: savings = 1x(ngridm) row vector of savings
		%				 shocks = (expn)x1 column vector of shocks
		%output: w1 = (expn)x(ngridm) matrix of all possible next period wealths
		w1=ones(size(shocks,1),1)*savings(1,:)*(1+me.r)+ ...
			 income(me,it+1,shocks(:,1))*ones(1,size(savings,2));
	end %budget
	function mw1=mbudget(me,it,savings,shocks,choice) 
		%derivative of wealth in t+1 w.r.t. savings
		%inputs and outputs as above
		mw1=repmat(1+me.r,size(shocks,1),size(savings,2));
	end %mbudget

	%Solver EGM
	function solve_dcegm(me)
		%solve the model with DC-EGM algorithm
		if me.lambda<eps & me.sigma>eps
			error(sprintf('Solution for the model with income shocks but without taste shocks is not supported!\nCan not use quadrature to calculate expectations of discontinuous functions.\nSet lambda > 0 or sigma = 0'))
		end
		me.policy=polyline;%delete previous solution
		me.value=polyline;%delete previous solution
		% first solve the retiree problem
		[quadp quadw]=quadpoints(me.expn,0,1);	%quadrature points
		quadstnorm=norminv(quadp,0,1);	%normally distributed
		savingsgrid=linspace(0,me.mmax,me.ngridm); %grid over savings (start with 0, ASCENDING)
    %main EGM loop
  	fprintf('t:');
    for it=me.Tbar:-1:1
    	fprintf(' %d',it);
	    if it==me.Tbar
        %terminal period
        for id=1:2 %1=staying, 2=leaving
	        me.policy(id,it)=polyline([0 me.mmax],[0 me.mmax],'Policy function in period T');
	        me.value(id,it)=polyline([0 me.mmax],[0 NaN],'Value function in period T');
	      	%vf(1)=0 is important, otherwise vf can not be computed in terminal period
	      end
	    else
        %not the terminal period
        for id=1:2 %for each decision %1=staying, 2=leaving
                wk1=me.budget(it,savingsgrid,quadstnorm*me.sigma, id);%next period wealth matrix
	            wk1=max(me.cfloor,wk1);%to insure minimum consumption and to avoid NaNs
	            vl1     =me.value_function(1,it+1,wk1(1:end));%next period value of staying, reshaped to row vector
	            vl1(2,:)=me.value_function(2,it+1,wk1(1:end));%next period value of leaving, row vector
	            pr1=me.chpr(vl1)*(id==1);%probability to stay tomorrow (conditional on being a worker)
	            cons1=me.policy(:,it+1).interpolate(wk1);%next period consumption by choice (2 rows reshaped matrix)
					mu1=pr1.*me.mutil(cons1(1,:))+(1-pr1).*me.mutil(cons1(2,:));%marginal utility by choice probs
	            mwk1=me.mbudget(it,savingsgrid,quadstnorm*me.sigma,id);%next period marginal wealth
	            rhs=quadw'*(reshape(mu1,size(wk1)).*mwk1);%RHS of Euler equation
	            cons0=me.imutil(me.df*rhs);%current period consumption
	            me.policy(id,it)=polyline(savingsgrid+cons0,cons0,sprintf('period %d, choice %d (1=staying)',it,id));
	          	%value function
	          	if id==1
	          		%stayer
	          		ev=quadw'*reshape(me.logsum(vl1),size(wk1));%expected value function for 1=staying
	          	else
	          		%emigrated
		            ev=quadw'*me.value_function(2,it+1,wk1);%expected value function for 2=leaving
          		end
          		me.value(id,it)=polyline(savingsgrid+cons0,me.util(cons0,id,it)+me.df*ev,sprintf('VF in period %d, choice %d (1=staying)',it,id));%uncleaned
          		%secondary envelope over bend over regions
          		if id==1 %only for stayers
	          		minx=min(me.value(id,it).x); %minimal x
	          		if me.value(id,it).x(1)<=minx
	          			% normal case: no bend back behind the first point
		          		[me.value(id,it) indxdel newdots]=me.value(id,it).secondary_envelope;
		          	else
		          		%special care needed in case the non-convex region coincides with credit constraint
		          		%this is when a bend back goes below the first point
		          		%value function to the left of the first point is analytical, so just add some points
		          		%there, so that polyline methods can be applied
		          		x1=linspace(minx,me.value(id,it).x(1),floor(me.ngridm/10)); %use .1 grid length
		          		x1(end)=[]; %skip the last point not to repeat it twice
		          		y1=me.util(x1,id)+me.df*ev(1); %analytical value function
		          		me.value(id,it)=me.value(id,it).inset(x1,y1,0);%insert all points in the front
		          		me.policy(id,it)=me.policy(id,it).inset(x1,x1,0);%keep the policy function on the same grid
		          		[me.value(id,it) indxdel newdots]=me.value(id,it).secondary_envelope;
	      	    	end
	      	    	me.value(id,it).label=sprintf('VF in period %d, choice %d (1=staying), cleaned',it,id);
			    			if numel(indxdel)>0
			    				%if any points were deleted in secondary envelope
			    				%analyze the policy function after the secondary envelope
			    				newpolicy=[];
			    				for i=1:numel(newdots.x)
			    					newpolicy(i,1)=newdots.x(i); %x for new point in policy
			    					j=find(me.policy(id,it).x<newdots.x(i)); %all point to the left of the new one
			    					j=j(~ismember(j,indxdel)); %not deleted
			    					j=max(j); %last such point
			    					newpolicy(i,2)=interp1(me.policy(id,it).x(j:j+1),me.policy(id,it).y(j:j+1),newdots.x(i),'linear');%interpolated from the left
			    					j=find(me.policy(id,it).x>newdots.x(i)); %all point to the right of the new one
			    					j=j(~ismember(j,indxdel)); %not deleted
			    					j=min(j); %first such point
			    					newpolicy(i,3)=interp1(me.policy(id,it).x(j-1:j),me.policy(id,it).y(j-1:j),newdots.x(i),'linear');%interpolated from the right
			    				end
			        		%remove inferior points from policy
			        		me.policy(id,it)=me.policy(id,it).thinout(indxdel);
			        		%add new points to policy twice to accurately describe discontinuities
			        		for i=1:size(newpolicy,1)
			        			j=find(me.policy(id,it).x>newpolicy(i,1),1,'first'); %first point past the insertion
			        			me.policy(id,it)=me.policy(id,it).inset(newpolicy(i,1)-1e3*eps,newpolicy(i,2),j-1);
			        			me.policy(id,it)=me.policy(id,it).inset(newpolicy(i,1),newpolicy(i,3),j);
			        		end
			    	    end
			    		end
	            %add special first point to create credit constrained region, save ev of saving zero
	      	    me.policy(id,it)=me.policy(id,it).inset(0,0,0);%connect the dots for the credit constrained region
	            me.value(id,it)=me.value(id,it).inset(0,ev(1),0);%last zero to make first point
	      end
	    end % if(terminal period)
    end %it
    fprintf('\n');
	end

    function res=solve_emigrated(me)
		%solve the emigrant problem with EGM algorithm
		%using res as output polyline
		res.policy=polyline;%delete previous solution
		res.value=polyline;%delete previous solution
		[quadp quadw]=quadpoints(me.expn,0,1);	%quadrature points
		quadstnorm=norminv(quadp,0,1);	%normally distributed
		savingsgrid=linspace(0,me.mmax,me.ngridm); %grid over savings (should start with 0)
		function vfres=vfcalc(it,x)
			vfres=nan(size(x)); %output of the same size as x
			mask=x<res.value(it).x(2); %all points in credit constrained region
			mask=mask | it==me.Tbar; %in the terminal period all points are in the constrained region
			vfres(mask)=me.util(x(mask),2)+me.df*res.value(it).y(1); %the first value in me.value(choice,it) is EV from zero savings!
			vfres(~mask)=res.value(it).interpolate(x(~mask));
		end
    %main EGM loop
    for it=me.Tbar:-1:1
	    if it==me.Tbar
        %terminal period
        res.policy(it)=polyline([0 me.mmax],[0 me.mmax],'Policy function in period T');
        res.value(it)=polyline([0 me.mmax],[0 NaN],'Value function in period T');
      	%vf(1)=0 is important, otherwise vf can not be computed in terminal period
	    else
        %not the terminal period
        wk1=me.budget(it,savingsgrid,quadstnorm*me.sigma,0);%next period wealth matrix
        cons1=res.policy(it+1).interpolate(wk1);%next period consumption
        cons1=max(cons1,me.cfloor);%safety net consumption floor
        mwk1=me.mbudget(it,savingsgrid,quadstnorm*me.sigma,0);%next period marginal wealth
        rhs=quadw'*(me.mutil(cons1).*mwk1);%RHS of Euler equation
        cons0=me.imutil(me.df*rhs);%current period consumption
        res.policy(it)=polyline(savingsgrid+cons0,cons0,sprintf('Policy in periond %d',it));
      	res.policy(it)=res.policy(it).inset(0,0,0);%connect the dots
        ev=quadw'*vfcalc(it+1,wk1);%expected value function
        res.value(it)=polyline(savingsgrid+cons0,me.util(cons0,0)+me.df*ev,sprintf('Value function in periond %d',it));
        %add special first point: ev of saving zero (unless already there)
        res.value(it)=res.value(it).inset(0,ev(1),0);%last zero to make first point
	    end % if(terminal period)
    end %it
	end


	%Solver VFI
	function solve_vfi(me)
		%solve the model with Value functions
		error 'Not implemented: homework'
	end

	%Solver Euler
	function solve_euler(me)
		%solve the model by solving Euler euqations
		error 'Not implemented: homework'
	end

	%Simulator
	function sim(me,seed)
		%simulate from the model
		%input: init = initial wealth
		%			  seed = seed for random number generator 
		%							(to run identical or varying simulations)
		if me.policy(1).len<1
			error 'The model should be solved first'
		end
		%fix the stream of random numbers
		if nargin<2
			rng(7134,'twister');
		else
			rng(seed,'twister');
		end
		%allocate
		me.sims=struct('wealth0',nan(me.nsims,me.Tbar), ...
									 'wealth1',nan(me.nsims,me.Tbar), ...
									 'consumption',nan(me.nsims,me.Tbar), ...
									 'shock',nan(me.nsims,me.Tbar), ...
									 'income',nan(me.nsims,me.Tbar), ...
									 'decision',nan(me.nsims,me.Tbar), ...
									 'prob_stay',nan(me.nsims,me.Tbar), ...
									 'migration_age',nan(me.nsims,1));
		%simulate choices and states
    for it=1:me.Tbar
    	if it==1
    		%draw initial wealth form uniform distribution on given interval
				me.sims.wealth0(:,it)=me.init(1)+rand(me.nsims,1)*(me.init(2)-me.init(1));
				me.sims.shock(:,it)=nan(me.nsims,1);
				me.sims.income(:,it)=nan(me.nsims,1);
				me.sims.decision(:,it)=ones(me.nsims,1); %everyone starts as a stayer
		else
				%choice state 
				me.sims.decision(choice,it)=1; %those staying remain stayers
				me.sims.decision(~choice,it)=2; %others leave
				me.sims.migration_age(me.sims.decision(:,it-1)==1 & me.sims.decision(:,it)==2)=it; %record retirement age
				%simulate wage shocks
				me.sims.shock(choice,it)=norminv(rand(sum(choice),1),0,1)*me.sigma;%normal scaled
                me.sims.shock(~choice,it)=norminv(rand(sum(~choice),1),0,1)*me.sigma;
				%compute wealth of all as emigrants first, then correct for
                %stayers
				%me.sims.wealth0(:,it)=diag(me.budget(it-1,me.sims.wealth1(:,it-1)',zeros(me.nsims,1),1));%use zero shock as it does not enter budget for retirees %%%%%%%%%%%%%%
                me.sims.wealth0(choice,it)=diag(me.budget(it-1,me.sims.wealth1(choice,it-1)',me.sims.shock(choice,it),1));%match savings and shocks one-to-one
                me.sims.wealth0(~choice,it)=diag(me.budget(it-1,me.sims.wealth1(~choice,it-1)',me.sims.shock(~choice,it),2));
                %me.sims.income(:,it)=0;
				me.sims.income(choice,it)=me.income(it,me.sims.shock(choice,it));
                me.sims.income(~choice,it)=me.income(it,me.sims.shock(~choice,it));
        end
	  	%choice probability to remain stayer
        vl1     =me.value_function(1,it,me.sims.wealth0(:,it)');
        vl1(2,:)=me.value_function(2,it,me.sims.wealth0(:,it)');
        me.sims.prob_stay(:,it)=(me.sims.decision(:,it)==1).*me.chpr(vl1)';%probability to stay (conditional on being a stayer)
        %simulate retirement choice
        choice=me.sims.prob_stay(:,it)>rand(me.nsims,1);
	  	%choice of consumption
	  	me.sims.consumption(choice,it)=me.policy(1,it).interpolate(me.sims.wealth0(choice,it));
	  	me.sims.consumption(~choice,it)=me.policy(2,it).interpolate(me.sims.wealth0(~choice,it));
	  	%end of period wealth
		me.sims.wealth1(:,it)=me.sims.wealth0(:,it)-me.sims.consumption(:,it);
    end
    end

	%Plotting
	function ax=plot(me,what2plot)
		%plot computed entities: policy, value, "sims smth"
		if nargin<2
			what2plot='solution';
		end
		what2plot=strsplit0(lower(what2plot),' ');%explode by space, convert to cell array
		if sum(ismember(what2plot,{'policy','solution','pol'}))>0
			%plotting the policy functions
			me.make_plot(me.policy,sprintf('%s: %s',me.label,'optimal consumption rules'));
		elseif sum(ismember(what2plot,{'value','value_function','val','valfunc','vf'}))>0
			%replace the analytical region with polylines
			k=100;%points to be added in analytical region
			[~,data]=me.value.chop(1); %through first points from all value functions
			for it=1:me.Tbar
				for id=1:2
					pt=exp(linspace(log(.01),log(me.value(id,it).x(2)),k));%log-grid
					data(id,it)=data(id,it).grow(polyline(pt,me.value_function(id,it,pt)),true); %true for adding in the front
					data(id,it).label=sprintf('t=%d id=%d',it,id);
				end
			end
			ax=me.make_plot(data,sprintf('%s: %s',me.label,'value functions'));
			ymax=max([data.y]);
			ymax=(floor(ymax/5)+1)*5; %round up to 5
			set(ax(1),'YLim',[-ymax ymax]);
			set(ax(2),'YLim',[-ymax ymax]);
		elseif sum(ismember(what2plot,{'prob','prob_stay','p','pr','pw','ch','ccp'}))>0
			%probability of staying (from values)
			%replace the analytical region with polylines
			k=100;%points to be added in analytical region
			[~,data]=me.value.chop(1); %through first points from all value functions
			for it=1:me.Tbar
				for id=1:2
					pt=exp(linspace(log(.01),log(me.value(id,it).x(2)),k));%log-grid
					data(id,it)=data(id,it).grow(polyline(pt,me.value_function(id,it,pt)),true); %true for adding in the front
				end
			end
			%compute choice probabilities and plot
			position=get(0,'ScreenSize');
			position([1 2])=0;position=position-[0 0 0 100];
			fig1 = figure('Color','white','Position',position);
			ax = axes('Parent',fig1,'FontSize',14);
			for it=1:me.Tbar
				%use data for stayers, and interpolate for emigrant
				xx=data(1,it).x;
				yy=me.chpr([data(1,it).y;data(2,it).x]);
				plot(xx,yy,'Parent',ax,'Marker','none','LineWidth',1,'DisplayName',sprintf('it=%d',it));
				hold(ax,'all');
			end
			set(ax,'XLim',[0 me.mmax],'YGrid','on','XGrid','on','YLim',[0 1]);
			box(ax,'on');
			xlabel(ax,'Wealth','Interpreter','None','FontSize',14);
			ylabel(ax,'Probability to stay by age','Interpreter','None','FontSize',14);
			drawnow
		elseif sum(ismember(what2plot,{'simulations','simulation','sim','sims'}))>0
			if numel(me.sims)==0
				error 'Nothing to plot'
			end
			flds=fields(me.sims)';
			map=ismember(flds,what2plot);
			map=map | sum(map)==0; %plot everything if nothing is chosen
			for fld=flds(map)
				if numel(me.sims.(fld{1}))<1
					error 'Nothing to plot'
				end
				if size(me.sims.(fld{1}),2)<2
					if size(me.sims.(fld{1}),1)>1
						%histogram if multiple sims
						position=get(0,'ScreenSize');
						position([1 2])=0;position=position-[0 0 0 100];
						fig1 = figure('Color','white','Position',position);
						ax = axes('Parent',fig1,'FontSize',14);
						histogram(me.sims.(fld{1}),'Parent',ax);
						title(ax,sprintf('%s: migration age distribution',me.label),'Interpreter','None','FontSize',14);
					end
				else
					%lifecycle plots
					position=get(0,'ScreenSize');
					position([1 2])=0;position=position-[0 0 0 100];
					fig1 = figure('Color','white','Position',position);
					ax = axes('Parent',fig1,'FontSize',14);
					h=plot([1:me.Tbar]',me.sims.(fld{1})');
					set(h,'Color','k','LineWidth',.5);
					hold(ax,'all');
					set(ax,'XLim',[1 me.Tbar],'YGrid','of','XGrid','of');
					%add migration markers
					xx=me.sims.migration_age-1;
					xx(isnan(xx))=me.Tbar;
					yy=diag(me.sims.(fld{1})(:,xx));
					yy(isnan(me.sims.migration_age))=NaN;
					stem(xx,yy,'MarkerSize',4,'Marker','square', ...
						'LineStyle','--','Color','k','DisplayName','Last year of stay');
					%prevent Y scale from getting too detailed
					corr=@(x,tol) [mean(x)-max(max(x)-min(x),tol)/2,mean(x)+max(max(x)-min(x),tol)/2];
					set(ax,'YLim',corr(get(ax,'YLim'),1e-8));
					box(ax,'off');
					xlabel(ax,'Age','Interpreter','None','FontSize',14);
					title(ax,sprintf('%s: simulated %s',me.label,fld{1}),'Interpreter','None','FontSize',14);
				end
			end
		else
			error 'Didn''t understand what to plot..'
		end
	end

	%Calculator of value functions that uses the analytical part in credit constrained region
	function res=value_function(me,choice,it,x)
		%interpolates value function at period t=it using analytical part
		%choice =1 for stayers, =2 for emigrants
		if me.value(choice,it).len<2
			error(sprintf('Can not compute value function at period %d because it only has %d points',it,me.value(choice,it).len))
		end
		res=nan(size(x)); %output of the same size as x
		mask=x<me.value(choice,it).x(2); %all points in credit constrained region
		mask=mask | it==me.Tbar; %in the terminal period all points are in the constrained region
		res(mask)=me.util(x(mask),choice,it)+me.df*me.value(choice,it).y(1); %the first value in me.value(choice,it) is EV from zero savings!
		res(~mask)=me.value(choice,it).interpolate(x(~mask));
	end

	%Logsum and choice probability calculators, assume numel(x)=2
	function res=logsum(me,x)
		%logsum by columns
		mx=max(x,[],1);
		mxx=x-repmat(mx,size(x,1),1);
		res=mx+me.lambda*log(sum(exp(mxx/me.lambda),1));
	end
	function res=chpr(me,x)
		%choice probability of the first row in multirow matrix
		mx=max(x,[],1);
		mxx=x-repmat(mx,size(x,1),1);
		res=exp(mxx(1,:)/me.lambda)./sum(exp(mxx/me.lambda),1);
	end

end %methods


methods (Access=private)
	function ax=make_plot(me,datalines,titlestr)
			if datalines(1).len==0
				error 'Nothing to plot'
			end
			position=get(0,'ScreenSize');
			position([1 2])=0;position=position-[0 0 0 100];
			fig1 = figure('Color','white','Position',position);
			ax(1) = subplot(1,2,1,'Parent',fig1,'FontSize',14);
			ax(2) = subplot(1,2,2,'Parent',fig1,'FontSize',14);
			datalines(1,:).plot(ax(1),'Marker','none','LineWidth',1);
			datalines(2,:).plot(ax(2),'Marker','none','LineWidth',1);
			for j=1:2
				hold(ax(j),'all');
				set(ax(j),'XLim',[0 me.mmax],'YGrid','on','XGrid','on');
				box(ax(j),'on');
				xlabel(ax(j),'Wealth','Interpreter','None','FontSize',14);
			end
			title(ax(1),'Stayer','Interpreter','None','FontSize',12);
			title(ax(2),'Emigrated','Interpreter','None','FontSize',12);
			mtit(fig1,titlestr,'Interpreter','None','FontSize',14,'yoff',0.05);
			drawnow
	end
end
end %of classdef

