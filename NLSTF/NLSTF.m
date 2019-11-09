
function HRHSI = NLSTF(HSI,MSI,F,downsampling_scale,para,par)
patchsize=8;
overlap=4;
 bparams.block_sz = [patchsize, patchsize];
 bparams.overlap_sz=[overlap overlap];
M=size(MSI,1);
N=size(MSI,2);
L=size(HSI,3);


num1=(M-patchsize)/(patchsize-overlap)+1;
num2=(N-patchsize)/(patchsize-overlap)+1;

 bparams.block_num=[num1 num2]; 
Z=zeros(patchsize,patchsize,L, bparams.block_num(1)* bparams.block_num(2));

  
 predenoised_blocks = ExtractBlocks(MSI, bparams); 
 Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);
 if para.K==1
     a=ones(num1*num2,1);
 else
%      gggggggg.careful=1;
   [a ]=fkmeans(Y2,para.K);
 end
  
 for mn=1:max(a)

 gg=find(a==mn);
 
 b=zeros(L,4*length(gg));
dddd=(downsampling_scale/4);
kkkk=1;
for i=1:length(gg)
     n=ceil(gg(i)/num2);
    m=gg(i)-num1*(n-1);
    
%     h1=ceil(m/dddd);
%     h2=ceil(n/dddd);
%     h3=min(h1+1,M/downsampling_scale);
%     h4=min(h2+1,N/downsampling_scale);
%       kk= HSI(h1:h3, h2:h4,:);
% kk1 = hyperConvert2D(kk);
%  pp=size(kk1,2);
%    b(:,kkkk:kkkk+pp-1)=kk1; 
%     kkkk=kkkk+pp;
    
    
    if (mod(m,dddd)~=0)&&(mod(n,dddd)~=0)

    kk= HSI(ceil(m/dddd), ceil(n/dddd),:);
    kk=kk(:);
    b(:,kkkk)=kk;
    kkkk=kkkk+1;
    end
    
     if (mod(m,dddd)==0)&&(mod(n,dddd)~=0)
       kk= HSI(ceil(m/dddd), ceil(n/dddd),:);
    kk=kk(:);
   b(:,kkkk)=kk;
    kkkk=kkkk+1;

     kk= HSI(ceil(m/dddd)+1, ceil(n/dddd),:);
       kk=kk(:);
 b(:,kkkk)=kk;
    kkkk=kkkk+1;
     
     end
      if (mod(m,dddd)~=0)&&(mod(n,dddd)==0)  
           kk= HSI(ceil(m/dddd), ceil(n/dddd),:);
    kk=kk(:);
    b(:,kkkk)=kk;
    kkkk=kkkk+1;
   
      kk= HSI(ceil(m/dddd), ceil(n/dddd)+1,:);
    kk=kk(:);
    b(:,kkkk)=kk;
    kkkk=kkkk+1;
     
      end
      if (mod(m,dddd)==0)&&(mod(n,dddd)==0) 
           kk= HSI(ceil(m/dddd), ceil(n/dddd),:);
    kk=kk(:);
      b(:,kkkk)=kk;
    kkkk=kkkk+1;
     kk= HSI(ceil(m/dddd)+1, ceil(n/dddd),:);
    kk=kk(:);
    b(:,kkkk)=kk;
    kkkk=kkkk+1;
     kk= HSI(ceil(m/dddd), ceil(n/dddd)+1,:);
    kk=kk(:);
     b(:,kkkk)=kk;
    kkkk=kkkk+1;
    kk= HSI(ceil(m/dddd)+1, ceil(n/dddd)+1,:);
    kk=kk(:);
    b(:,kkkk)=kk;
    kkkk=kkkk+1;
      end
end
 b=b(:,1:kkkk-1);
  zz=unique(b','rows');
  zz=zz';
[~, nn]=size(zz);
 par.K=min(para.S,nn);
  par.lambda=para.lambda3;
%    phi  =    Nonnegative_DL( zz, par );
     phi=vca(zz,par.K);
%    phi= VCA(zz,'Endmembers',par.K,'SNR',0,'verbose','off');


% %  phi=phi(:,1:min(par.K,size(phi,2)));
% max_vol = 0;
%     vol = zeros(1, 10);
%     for idx_VCA = 1:10
% %          phi1= VCA(zz,'Endmembers',par.K,'SNR',0,'verbose','off');
% phi11=vca(zz,par.K);
%         vol(idx_VCA) = abs(det(phi11'*phi11));
%         if vol(idx_VCA) > max_vol
%             phi = phi11;
%             max_vol = vol(idx_VCA);
%         end   
%     end
         ph=F*phi;

      
       %ÑµÁ·×Öµäphi1
          BB=predenoised_blocks(:,:,:,gg);
      cc=Unfold(BB,size(BB),1);
   cc=unique(cc','rows');
   cc=cc';
         
 D = cc;
   [~ ,nn]=size(cc);
    par.lambda=para.lambda1;
 par.K=min(para.W,nn);
phi1  =    Nonnegative_DL( D, par );
       
          %ÑµÁ·×Öµäphi2
           
      cc=Unfold(BB,size(BB),2);
       cc=unique(cc','rows');
   cc=cc';
 D = cc;
   [~,nn]=size(cc);
        par.lambda=para.lambda2;
 par.K=min(para.H,nn);
phi2  =    Nonnegative_DL( D, par );

  D1=phi1;
  D2=phi2;
  D3=ph;

D9{1}=D1;
D9{2}=D2;
D9{3}=D3;
  x   = sparse_tucker( D9, BB, para.lambda );
%  x   = l2_tucker( D9, BB, para.lambda );
if ismatrix(x)
    d=ttm(tensor(x),{D1,D2},[1,2]);
   d=permute(double(d),[3 1 2]);
   d=ttm(tensor(d),{phi},1);
     d=permute(double(d),[2 3 1]);
else
d=ttm(tensor(x),{D1,D2,phi},[1,2,3]);
end
d=double(d);
 Z(:,:,:,gg)=d;
end

EZ = JointBlocks(Z, bparams);

G=create_G([1 M/downsampling_scale 1 N/downsampling_scale], downsampling_scale);
b11=alternating_back_projection(hyperConvert2D(EZ),hyperConvert2D(HSI),hyperConvert2D(MSI),F,G,par);
HRHSI=hyperConvert3D(b11,M,N);




 
