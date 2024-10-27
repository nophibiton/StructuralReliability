function [c_sys,sys_type,n_cut_sets]=sys_event(sys_def,C)

n=length(sys_def{1});
[row,col]=size(C);

if n==1 && sys_def{1}<0;
    sys_type='series';
    c_sys=1-prod(ones(row,col)-C,2);
    n_cut_sets=[];   
elseif n==1 && sys_def{1}>0;
    sys_type='parallel';
    c_sys=prod(C,2);
    n_cut_sets=[];    
else
    sys_type='general';
    if strcmpi(sys_def{2},'link')==1; sys_def{1}=-sys_def{1}; end;
    sys_def{1}=[0 sys_def{1} 0];
    sys_nz=find(sys_def{1}~=0); sys_z=find(sys_def{1}==0);
    int_1=sys_z-[0 sys_z(1:end-1)];
    size_cut_sets=int_1(int_1>1)-1; n_cut_sets=length(size_cut_sets);    
    for i=1:n_cut_sets;
        c_cut_set=ones(row,1);
        for j=1:size_cut_sets(i);
            comp=sys_def{1}(sys_nz(sum(size_cut_sets(1:i-1))+j));
            if comp<0;
                c_cut_set=c_cut_set.*(ones(row,1)-C(:,abs(comp)));
            else
                c_cut_set=c_cut_set.*C(:,comp);
            end
        end
        c_cut_sets(:,i)=c_cut_set;
    end
    c_sys=ones(row,1)-prod(ones(row,n_cut_sets)-c_cut_sets,2);    
end

end
