Num_rep=10;
Dim_array=[5,7,9,10,100,200];
Vtx_array=[50,100,200,300,500,1000];

frac_array=[0,20/100,50/100];
time_ta=zeros(length(frac_array),length(Dim_array));
time_scale_ta=zeros(length(frac_array),length(Dim_array));
time_asfw=zeros(length(frac_array),length(Dim_array));
time_fw=zeros(length(frac_array),length(Dim_array));
problem_size=zeros(1,length(frac_array));
for ii=1:length(Dim_array)
    Num_dim=Dim_array(ii);
    

    for jj=1:length(frac_array)
        frac_rat=frac_array(jj);
        Num_vtx=ceil(Vtx_array(ii)*(1-frac_rat))
        Num_pts=ceil(Vtx_array(ii)*(frac_rat))
        for kk=1:Num_rep
            ii
            jj
            kk
            vtx_A=Random_pts(Num_dim,Num_vtx,'unit ball');
            if frac_rat~=0
                inhulldata=Random_cvx(vtx_A,Num_pts,'dir');
                matA=[inhulldata,vtx_A];
                [m,n]=size(matA);
                normal_val_index=(Num_pts+1):n;

                rndindx=randperm(n,n);
                [or_val or_ind]=sort(rndindx);
                vertices_index=or_ind(normal_val_index);
                matA=matA(:,rndindx);
            else
                matA=vtx_A;
                [m,n]=size(matA);
            end
            disp('ta')
            tic;
            [index_this_avt,iter_ta_val]=AVTA_anti(matA',0.0005);
            ta_end=toc;
            time_ta(jj,ii)=time_ta(jj,ii)+ta_end;
            disp('end ta')
            disp('spta')
            tic;
            [index_this_spavt,iter_spta]=AVTA_SP(matA',0.0005);
            scale_ta_end=toc;
            time_scale_ta(jj,ii)=time_scale_ta(jj,ii)+scale_ta_end;
            disp('end spta')
            disp('FW')
            tic;
            [index_this_fw,iter_fw]=AVTA_GT(matA',0.0005);
            fw_end=toc;
            time_fw(jj,ii)=time_fw(jj,ii)+fw_end;
            disp('end FW')
            disp('ASFW')
            tic;
            [index_this_asfw,iter_asfw]=AVTA_ASFW(matA',0.0005);
            asfw_end=toc;
            time_asfw(jj,ii)=time_asfw(jj,ii)+asfw_end;
            disp('end FW')
         end
        
    end
    
    
end
        






names=strings(6,1);
for i=1:length((time_ta(1,:)))
    this_size=[num2str(Dim_array(i)), 'x', num2str(Vtx_array(i))];
    names(i)=this_size;
end
        
for i=1:length((time_ta(:,1)))

    this_title=['vertice generated on Unit Sphere ', num2str(frac_array(i)*100),'% redundant point'];
    figure(i);
    hold on
    title(this_title);

    plot(log(time_ta(i,:)),'DisplayName','TA','LineWidth',1.5)
    plot(log(time_scale_ta(i,:)),'DisplayName','SPH','LineWidth',1.5)  
    plot(log(time_fw(i,:)),'DisplayName','FW','LineWidth',1.5)  
    plot(log(time_asfw(i,:)),'DisplayName','ASFW','LineWidth',1.5)  
    legend('show','Location','northwest')%,'Orientation','horizontal')
    set(gca,'xtick',[1:length((time_ta(i,:)))],'xticklabel',names)
    xlabel ("Problem Size");
    ylabel ("Running time (secs in log scale)");
    save_title= strrep(this_title,' ','_');
    saveas(gcf,['allvertices_',save_title,'.png'])
    hold off;
end
