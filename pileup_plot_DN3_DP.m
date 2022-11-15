clc
clf
clear


ss = "DN3";
% ss = "DP";

lgene0 = 100000;
dx = 10000;
lg0 = 10;
lg1 = 40;
r = 0;



B = [197195432 181748087 159599783 155630120 152537259 149517037 152524553 131738871 124076172 129993255 121843856 121257530 120284312 125194864 103494974 98319150 95272651 90772031 61342430 166650296 91744698];
chrm = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"};


Hgpm = zeros(1,2*lg1+lg0,2*lg1+lg0);
Hgmm = zeros(1,2*lg1+lg0,2*lg1+lg0);
Hgp0m = zeros(1,2*lg1+lg0,2*lg1+lg0);
Hgm0m = zeros(1,2*lg1+lg0,2*lg1+lg0);
Hgtot = zeros(1,2*lg1+lg0,2*lg1+lg0);
ngp = 0;
ngm = 0;
genedata = zeros(1,1);

for ni=1:20
    ni
    chr = chrm{ni};
    ll = fix(B(ni)/dx);
    H = zeros(ll+1,ll+1);
    
    A = load(sprintf('data/HiC_%s_10kb/%s.txt',ss,chr));

    m = length(A);

    for i=1:m
        a = A(i,1);
        b = A(i,2);
        c = A(i,4);
        H(a,b) = c;
        H(b,a) = c;
    end  
    
    H0 = zeros(ll+1,ll+1);
    for i=0:ll
        s = zeros(1,1);
        k = 0;
        for j=1:ll-i
            if(H(j,j+i)>0)
                k = k + 1;
                s(k) =  H(j,j+i);
            end
        end
        sm = median(s);
        for j=1:ll-i
            H0(j,j+i) = sm;
            H0(j+i,j) = sm;
        end
    end
            
    RNAPII = load(sprintf('data/ChIP_PolII_%s/PolII_%s.txt',ss,chrm{ni}));
    genes0 = load(sprintf('data/gene_list_RNAseq_%s/gene_list_RNAseq_%s.txt',ss,chrm{ni}));
    k = 0;
    genes = zeros(1,1);
    for i= 1:length(genes0)-1
        if (genes0(i,1)>lgene0)
            k = k +1;
            for j=1:5
                genes(k,j) = genes0(i,j);
            end
        end
    end
    

    gene2 = genes;
    
    l = fix((B(ni))/dx)+1;
    RNAPII = sortrows(RNAPII);
    lRNA = length(RNAPII);
    x = zeros(l,1);
    j=1;
    for i=1:l-1
        n = 0;
        while (RNAPII(j,1)<=i*dx && j<lRNA)
            n = n + RNAPII(j,2);
            j = j + 1;
        end
        x(i) = n;
    end
    
    m1 = size(gene2);
    for i=1:m1(1)
        xi = fix(gene2(i,2)/dx)+1;
        xf = fix(gene2(i,3)/dx)+1;
        if (gene2(i,4)==1)

            tss(i) = x(xi);
            tts(i) = x(xf);
            
            s = zeros(1,1);
            k = 0;
            for j=xi+(r+1):xf-1
                k = k + 1;
                s(k) =  x(j);
            end
            dgene(i) = median(s)/1;
            
            s = zeros(1,1);
            k = 0;
            for j=(2*xi-xf-1):xi-2
                if (j>0)
                    k = k + 1;
                    s(k) = x(j);
                end
            end
            taill(i) = median(s);

            pausei(i) = tss(i)/dgene(i);

               xi = fix(gene2(i,2)/dx)+1;
               xf = fix(gene2(i,3)/dx)+1;
               lgene = xf-xi+1;
                m1 = mod(lgene,lg0);
                m2 = fix(lgene/lg0);
                
                s = zeros(1,1);
                kk = 0;
                for j=1:lgene
                    for k=xi-lg1:xf-j-lg1
                        if (H(k+j,k)>0)
                            kk = kk + 1;
                            s(kk) = H(k+j,k)/H0(k+j,k);
                        end
                    end
                end
                domain_hic = median(s);
                
                s = zeros(1,1);
                kk = 0;
                for j=0:lgene
                    for k=xi:xf-j
                        if (H(k+j,k)>0)
                            kk = kk + 1;
                            s(kk) = H(k+j,k)/H0(k+j,k);
                        end
                    end
                end
                intragene_hic = median(s);
                
                if (m2>0)
                    l = lgene+2*lg1;
                    H1 = zeros(2,l,l);
                    lh = length(H);
                    for j=1:l
                        for k =1:l
                            if ((xi+j-lg1-1)>0 && (xi+j-lg1-1)<=lh)
                                if ((xi+k-lg1-1)>0 && (xi+k-lg1-1)<=lh)
                                    if(H(xi+j-lg1-1,xi+k-lg1-1)>0)
                                        H1(1,j,k) = H(xi+j-lg1-1,xi+k-lg1-1);
                                        H1(2,j,k) = H(xi+j-lg1-1,xi+k-lg1-1)/H0(xi+j-lg1-1,xi+k-lg1-1);
                                    end
                                end
                            end
                        end
                    end
                    sH = zeros(lgene,2);
                    H2 = zeros(2,2*lg1+m2*lg0,2*lg1+m2*lg0);

                    if (m1>0)
                        for j=1:lgene
                            s = 0;
                            for k=1:l
                                if (H1(2,lg1+j,k)>0)
                                    s = s + log2(H1(2,lg1+j,k));
                                end
                            end
                            sH(j,1) = s;
                            sH(j,2) = j;
                        end
                        sH1 = sortrows(sH);
                        rb = zeros(m1,1);
                        for j=1:m1
                            rb(j) = sH1(j,2);
                        end
                        rb = sort(rb);

                        nj = 0;
                        nk = 0;
                        for j=1:l
                            c = 1;
                            if (j>=rb(1)+lg1 && j<=rb(m1)+lg1)
                                for f=1:m1
                                    if (j==rb(f)+lg1)
                                        c = 0;
                                        break;
                                    end
                                end
                            end
                            if (c==1)
                                nj = nj + 1;
                                nk = 0;
                                for k=1:j
                                    c1 = 1;
                                    if (k>=rb(1)+lg1 && k<=rb(m1)+lg1)
                                        for f=1:m1
                                           if (k==rb(f)+lg1) 
                                               c1 = 0;
                                               break
                                           end
                                        end
                                    end
                                    if (c1 == 1)
                                        nk = nk + 1;
                                        H2(1,nj,nk) = H1(1,j,k);
                                        H2(1,nk,nj) = H1(1,j,k);
                                        H2(2,nj,nk) = H1(2,j,k);
                                        H2(2,nk,nj) = H1(2,j,k);
                                    end
                                end
                            end
                        end
                    else
                        H2 = H1;
                    end
                    H3 = zeros(2,2*lg1+lg0,2*lg1+lg0);
                    if (m2>1)
                        for j=lg1+1:lg1+lg0
                                for k=lg1+1:j
                                    s1 = zeros(1,1);
                                    s2 = zeros(1,1);
                                    ns = 0;
                                    for f=1:m2
                                        for h=1:m2
                                            if(H2(2,lg1+m2*(j-lg1-1)+f,lg1+m2*(k-lg1-1)+h)>0)
                                                ns = ns + 1;
                                                s1(ns) = H2(1,lg1+m2*(j-lg1-1)+f,lg1+m2*(k-lg1-1)+h);
                                                s2(ns) = H2(2,lg1+m2*(j-lg1-1)+f,lg1+m2*(k-lg1-1)+h);
                                            end
                                        end
                                    end
                                    H3(1,j,k) = median(s1);
                                    H3(1,k,j) = median(s1);
                                    H3(2,j,k) = median(s2);
                                    H3(2,k,j) = median(s2);
                                end
                            end
                            for j=1:lg1
                                for k=lg1+1:lg1+lg0
                                    s1 = zeros(1,1);
                                    s2 = zeros(1,1);
                                    ns = 0;
                                    for h=1:m2
                                        if(H2(2,j,lg1+m2*(k-lg1-1)+h)>0)
                                            ns = ns + 1;
                                            s1(ns) = H2(1,j,lg1+m2*(k-lg1-1)+h);
                                            s2(ns) = H2(2,j,lg1+m2*(k-lg1-1)+h);
                                        end
                                    end
                                    H3(1,j,k) = median(s1);
                                    H3(1,k,j) = median(s1);
                                    H3(2,j,k) = median(s2);
                                    H3(2,k,j) = median(s2);
                                end
                            end
                            for j=lg1+1:lg1+lg0
                                for k=1:lg1
                                    s1 = zeros(1,1);
                                    s2 = zeros(1,1);
                                    ns = 0;
                                    for f=1:m2
                                        if(H2(2,lg1+m2*(j-lg1-1)+f,k)>0)
                                            ns = ns + 1;
                                            s1(ns) = H2(1,lg1+m2*(j-lg1-1)+f,k);
                                            s2(ns) = H2(2,lg1+m2*(j-lg1-1)+f,k);
                                        end
                                    end
                                    H3(1,j,k) = median(s1);
                                    H3(1,k,j) = median(s1);
                                    H3(2,j,k) = median(s2);
                                    H3(2,k,j) = median(s2);
                                end
                            end
                            for j=1:lg1
                               for k=1:lg1
                                    H3(1,j,k) = H2(1,j,k); 
                                    H3(2,j,k) = H2(2,j,k);
                               end
                           end
                           for j=1:lg1
                               for k=1:lg1
                                    H3(1,j+lg1+lg0,k+lg1+lg0) = H2(1,j+lg1+m2*lg0,k+lg1+m2*lg0); 
                                    H3(2,j+lg1+lg0,k+lg1+lg0) = H2(2,j+lg1+m2*lg0,k+lg1+m2*lg0);
                               end
                           end
                           for j=1:lg1
                               for k=1:lg1
                                    H3(1,j,k+lg1+lg0) = H2(1,j,k+lg1+m2*lg0);  
                                    H3(1,k+lg1+lg0,j) = H2(1,j,k+lg1+m2*lg0); 
                                    H3(2,j,k+lg1+lg0) = H2(2,j,k+lg1+m2*lg0);  
                                    H3(2,k+lg1+lg0,j) = H2(2,j,k+lg1+m2*lg0); 
                               end
                           end
                            
                    else
                        H3 = H2;
                    end
                    ngp = ngp + 1;
                    for j=1:2*lg1+lg0
                        for k=1:2*lg1+lg0
                            Hgpm(ngp,j,k) =  H3(1,j,k);
                            Hgp0m(ngp,j,k) =  H3(2,j,k);
                        end
                    end
                    
                    
                    for j=1:2*lg1+lg0
                        for k=1:2*lg1+lg0
                            Hgtot(ngp+ngm,j,k) =  H3(2,j,k);
                        end
                    end
                    genedata(ngp+ngm,1) = tss(i);
                    genedata(ngp+ngm,2) = tts(i);
                    genedata(ngp+ngm,3) = dgene(i);
                    genedata(ngp+ngm,4) = gene2(i,1);
                    genedata(ngp+ngm,5) = gene2(i,5);
                    genedata(ngp+ngm,6) = gene2(i,4);
                    genedata(ngp+ngm,7) = domain_hic;
                    genedata(ngp+ngm,8) = intragene_hic;
                   
                else
                    m2 = 1;
                    lgene1 = fix(lgene/2);
                    lgene2 = lg0 -lgene+lgene1+1;
                    l = lgene+2*lg0;
                    l1 = lg0+lgene1;
                    l2 = lg0+lgene2;
                    H1 = zeros(2,3*lg0,3*lg0);
                    lh = length(H);
                    j1 = 0;
                    k1 = 0;
                    for j=1:3*lg0
                        if (j<=l1 || j>=l2)
                           j1 = j1 + 1; 
                           k1= 0;
                            for k =1:3*lg0
                                if (k<=l1 || k>=l2)
                                    k1 = k1 + 1;
                                    if ((xi+j1-lg0-1)>0 && (xi+j1-lg0-1)<=lh)
                                        if ((xi+k1-lg0-1)>0 && (xi+k1-lg0-1)<=lh)
                                            if(H(xi+j1-lg0-1,xi+k1-lg0-1)>0)
                                                H1(1,j,k) = H(xi+j1-lg0-1,xi+k1-lg0-1);
                                                H1(2,j,k) = H(xi+j1-lg0-1,xi+k1-lg0-1)/H0(xi+j1-lg0-1,xi+k1-lg0-1);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    sH = zeros(lgene,2);
                    H3 = H1;

                    for j=0:lg0-1              
                        s1 = zeros(1,1);
                        s2 = zeros(1,1);
                        ns = 0;            
                        for k=1:lg0-j
                            if (H3(2,lg0+k,lg0+k+j)>0)
                                ns = ns + 1;
                                s1(ns) =  H3(1,lg0+k,lg0+k+j);
                                s2(ns) =  H3(2,lg0+k,lg0+k+j);
                            end
                        end
                        if (ns>0)
                            ms1 = median(s1);
                            ms2 = median(s2);
                            for k=1:lg0-j
                                if (H3(2,lg0+k,lg0+k+j)<=0)
                                    H3(1,lg0+k,lg0+k+j)=ms1;
                                    H3(2,lg0+k,lg0+k+j)=ms2;
                                    H3(1,lg0+k+j,lg0+k)=ms1;
                                    H3(2,lg0+k+j,lg0+k)=ms2;
                                end
                            end
                        end    
                    end
                    ns = 0;
                        s1 = zeros(1,1);
                        s2 = zeros(1,1);
                        for j=1:lg0
                            for k=1:lg0
                                if (H3(2,lg0+j,lg0+k)>0)
                                    ns = ns + 1;
                                    s1(ns) = H3(1,lg0+j,lg0+k);
                                    s2(ns) = H3(2,lg0+j,lg0+k);
                                end
                            end
                        end
                        if (ns>0)
                            ms1 = median(s1);
                            ms2 = median(s2);
                            for j=1:lg0
                                for k=1:lg0
                                    if (H3(2,lg0+j,lg0+k)<=0)
                                        H3(1,lg0+j,lg0+k)=ms1;
                                        H3(2,lg0+j,lg0+k)=ms2;
                                    end
                                end
                            end
                        end
                        
                        for j=1:lg0
                            ns = 0;
                            s1 = zeros(1,1);
                            s2 = zeros(1,1);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)>0)
                                    ns = ns + 1;
                                    s1(ns) = H3(1,k,j);
                                    s2(ns) = H3(2,k,j);
                                end
                            end
                            ms1 = median(s1);
                            ms2 = median(s2);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)==0)
                                    sm1 = (ms1/2)*(0.5*rand()+1.5);
                                    sm2 = (ms2/2)*(0.5*rand()+1.5);
                                    H3(1,k,j)=sm1;
                                    H3(1,j,k)=sm1;
                                    H3(2,k,j)=sm2;
                                    H3(2,j,k)=sm2;
                                end
                            end
                        end
                        for j=2*lg0+1:3*lg0
                            ns = 0;
                            s1 = zeros(1,1);
                            s2 = zeros(1,1);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)>0)
                                    ns = ns + 1;
                                    s1(ns) = H3(1,k,j);
                                    s2(ns) = H3(2,k,j);
                                end
                            end
                            ms1 = median(s1);
                            ms2 = median(s2);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)==0)
                                    sm1 = (ms1/2)*(0.5*rand()+1.5);
                                    sm2 = (ms2/2)*(0.5*rand()+1.5);
                                    H3(1,k,j)=sm1;
                                    H3(1,j,k)=sm1;
                                    H3(2,k,j)=sm2;
                                    H3(2,j,k)=sm2;
                                end
                            end
                        end
                    ngp = ngp + 1;
                    for j=1:3*lg0
                        for k=1:3*lg0
                            Hgpm(ngp,j,k) =  H3(1,j,k);
                            Hgp0m(ngp,j,k) =  H3(2,j,k);
                        end
                    end
                    
                    
                end
            
            
        else

            tss(i) = x(xf);
            tts(i) = x(xi);
            
            s = zeros(1,1);
            k = 0;
            for j=xi+1:xf-(r+1)
                k = k + 1;
                s(k) =  x(j);
            end
            dgene(i) = median(s)/1;
            
            s = zeros(1,1);
            k = 0;
            for j=(xf+1):(2*xf-xi+1)
                if(j<=ll)
                    k = k + 1;
                    s(k) = x(j);
                end
            end
            taill(i) = median(s);

            pausei(i) = tss(i)/dgene(i);

                    xi = fix(gene2(i,2)/dx)+1;
                    xf = fix(gene2(i,3)/dx)+1;
                    lgene = xf-xi+1;
                    
                    s = zeros(1,1);
                    kk = 0;
                    for j=1:lgene
                        for k=xi-lgene:xf-j-lgene
                            if (H(k+j,k)>0)
                                kk = kk + 1;
                                s(kk) = H(k+j,k)/H0(k+j,k);
                            end
                        end
                    end
                    domain_hic = median(s);
                    
                    s = zeros(1,1);
                    kk = 0;
                    for j=0:lgene
                        for k=xi:xf-j
                            if (H(k+j,k)>0)
                                kk = kk + 1;
                                s(kk) = H(k+j,k)/H0(k+j,k);
                            end
                        end
                    end
                    intragene_hic = median(s);
                
                    m1 = mod(lgene,lg0);
                    m2 = fix(lgene/lg0);
                    if(m2>0)
                        l = lgene+2*lg1;
                    H1 = zeros(2,l,l);
                    lh = length(H);
                    for j=1:l
                        for k =1:l
                            if ((xi+j-lg1-1)>0 && (xi+j-lg1-1)<=lh)
                                if ((xi+k-lg1-1)>0 && (xi+k-lg1-1)<=lh)
                                    if(H(xi+j-lg1-1,xi+k-lg1-1)>0)
                                        H1(1,j,k) = H(xi+j-lg1-1,xi+k-lg1-1);
                                        H1(2,j,k) = H(xi+j-lg1-1,xi+k-lg1-1)/H0(xi+j-lg1-1,xi+k-lg1-1);
                                    end
                                end
                            end
                        end
                    end
                    sH = zeros(lgene,2);
                    H2 = zeros(2,2*lg1+m2*lg0,2*lg1+m2*lg0);

                    if (m1>0)
                        for j=1:lgene
                            s = 0;
                            for k=1:l
                                if (H1(2,lg1+j,k)>0)
                                    s = s + log2(H1(2,lg1+j,k));
                                end
                            end
                            sH(j,1) = s;
                            sH(j,2) = j;
                        end
                        sH1 = sortrows(sH);
                        rb = zeros(m1,1);
                        for j=1:m1
                            rb(j) = sH1(j,2);
                        end
                        rb = sort(rb);

                        nj = 0;
                        nk = 0;
                        for j=1:l
                            c = 1;
                            if (j>=rb(1)+lg1 && j<=rb(m1)+lg1)
                                for f=1:m1
                                    if (j==rb(f)+lg1)
                                        c = 0;
                                        break;
                                    end
                                end
                            end
                            if (c==1)
                                nj = nj + 1;
                                nk = 0;
                                for k=1:j
                                    c1 = 1;
                                    if (k>=rb(1)+lg1 && k<=rb(m1)+lg1)
                                        for f=1:m1
                                           if (k==rb(f)+lg1) 
                                               c1 = 0;
                                               break
                                           end
                                        end
                                    end
                                    if (c1 == 1)
                                        nk = nk + 1;
                                        H2(1,nj,nk) = H1(1,j,k);
                                        H2(1,nk,nj) = H1(1,j,k);
                                        H2(2,nj,nk) = H1(2,j,k);
                                        H2(2,nk,nj) = H1(2,j,k);
                                    end
                                end
                            end
                        end
                    else
                        H2 = H1;
                    end
                    H3 = zeros(2,2*lg1+lg0,2*lg1+lg0);
                        if (m2>1)
                            for j=lg1+1:lg1+lg0
                                for k=lg1+1:j
                                    s1 = zeros(1,1);
                                    s2 = zeros(1,1);
                                    ns = 0;
                                    for f=1:m2
                                        for h=1:m2
                                            if(H2(2,lg1+m2*(j-lg1-1)+f,lg1+m2*(k-lg1-1)+h)>0)
                                                ns = ns + 1;
                                                s1(ns) = H2(1,lg1+m2*(j-lg1-1)+f,lg1+m2*(k-lg1-1)+h);
                                                s2(ns) = H2(2,lg1+m2*(j-lg1-1)+f,lg1+m2*(k-lg1-1)+h);
                                            end
                                        end
                                    end
                                    H3(1,j,k) = median(s1);
                                    H3(1,k,j) = median(s1);
                                    H3(2,j,k) = median(s2);
                                    H3(2,k,j) = median(s2);
                                end
                            end
                            for j=1:lg1
                                for k=lg1+1:lg1+lg0
                                    s1 = zeros(1,1);
                                    s2 = zeros(1,1);
                                    ns = 0;
                                        for h=1:m2
                                            if(H2(2,j,lg1+m2*(k-lg1-1)+h)>0)
                                                ns = ns + 1;
                                                s1(ns) = H2(1,j,lg1+m2*(k-lg1-1)+h);
                                                s2(ns) = H2(2,j,lg1+m2*(k-lg1-1)+h);
                                            end
                                        end
                                    H3(1,j,k) = median(s1);
                                    H3(1,k,j) = median(s1);
                                    H3(2,j,k) = median(s2);
                                    H3(2,k,j) = median(s2);
                                end
                            end
                            for j=lg1+1:lg1+lg0
                                for k=1:lg1
                                    s1 = zeros(1,1);
                                    s2 = zeros(1,1);
                                    ns = 0;
                                    for f=1:m2
                                            if(H2(2,lg1+m2*(j-lg1-1)+f,k)>0)
                                                ns = ns + 1;
                                                s1(ns) = H2(1,lg1+m2*(j-lg1-1)+f,k);
                                                s2(ns) = H2(2,lg1+m2*(j-lg1-1)+f,k);
                                            end
                                    end
                                    H3(1,j,k) = median(s1);
                                    H3(1,k,j) = median(s1);
                                    H3(2,j,k) = median(s2);
                                    H3(2,k,j) = median(s2);
                                end
                            end
                            for j=1:lg1
                               for k=1:lg1
                                        H3(1,j,k) = H2(1,j,k); 
                                        H3(2,j,k) = H2(2,j,k);
                               end
                           end
                           for j=1:lg1
                               for k=1:lg1
                                        H3(1,j+lg1+lg0,k+lg1+lg0) = H2(1,j+lg1+m2*lg0,k+lg1+m2*lg0); 
                                        H3(2,j+lg1+lg0,k+lg1+lg0) = H2(2,j+lg1+m2*lg0,k+lg1+m2*lg0);
                               end
                           end
                           for j=1:lg1
                               for k=1:lg1
                                        H3(1,j,k+lg1+lg0) = H2(1,j,k+lg1+m2*lg0);  
                                        H3(1,k+lg1+lg0,j) = H2(1,j,k+lg1+m2*lg0); 
                                        H3(2,j,k+lg1+lg0) = H2(2,j,k+lg1+m2*lg0);  
                                        H3(2,k+lg1+lg0,j) = H2(2,j,k+lg1+m2*lg0); 
                               end
                            end
                        else
                            H3 = H2;
                        end
                        ngm = ngm + 1;
                        for j=1:2*lg1+lg0
                            for k=1:2*lg1+lg0
                                Hgmm(ngm,j,k) =  H3(1,j,k);
                                Hgm0m(ngm,j,k) =  H3(2,j,k);
                            end
                        end


                        for j=1:2*lg1+lg0
                            for k=1:2*lg1+lg0
                                Hgtot(ngp+ngm,j,k) =  H3(2,j,k);
                            end
                        end
                            genedata(ngp+ngm,1) = tss(i);
                            genedata(ngp+ngm,2) = tts(i);
                            genedata(ngp+ngm,3) = dgene(i);
                            genedata(ngp+ngm,4) = gene2(i,1);
                            genedata(ngp+ngm,5) = gene2(i,5);
                            genedata(ngp+ngm,6) = gene2(i,4);
                            genedata(ngp+ngm,7) = domain_hic;
                            genedata(ngp+ngm,8) = intragene_hic;

                        else
                        m2 = 1;
                        lgene1 = fix(lgene/2);
                        lgene2 = lg0 -lgene+lgene1+1;
                        l = lgene+2*lg0;
                        l1 = lg0+lgene1;
                        l2 = lg0+lgene2;
                        H1 = zeros(2,3*lg0,3*lg0);
                        lh = length(H);
                        j1 = 0;
                        k1 = 0;
                        for j=1:3*lg0
                            if (j<=l1 || j>=l2)
                               j1 = j1 + 1; 
                               k1 = 0;
                                for k =1:3*lg0
                                    if (k<=l1 || k>=l2)
                                        k1 = k1 + 1;
                                        if ((xi+j1-lg0-1)>0 && (xi+j1-lg0-1)<=lh)
                                            if ((xi+k1-lg0-1)>0 && (xi+k1-lg0-1)<=lh)
                                                if(H(xi+j1-lg0-1,xi+k1-lg0-1)>0)
                                                    H1(1,j,k) = H(xi+j1-lg0-1,xi+k1-lg0-1);
                                                    H1(2,j,k) = H(xi+j1-lg0-1,xi+k1-lg0-1)/H0(xi+j1-lg0-1,xi+k1-lg0-1);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        sH = zeros(lgene,2);
                        H3 = H1;
                        for j=0:lg0-1              
                            s1 = zeros(1,1);
                            s2 = zeros(1,1);
                            ns = 0;            
                            for k=1:lg0-j
                                if (H3(1,lg0+k,lg0+k+j)>0)
                                    ns = ns + 1;
                                    s1(ns) =  H3(1,lg0+k,lg0+k+j);
                                    s2(ns) =  H3(2,lg0+k,lg0+k+j);
                                end
                            end
                            if (ns>0)
                                ms1 = median(s1);
                                ms2 = median(s2);
                                for k=1:lg0-j
                                    if (H3(1,lg0+k,lg0+k+j)<=0)
                                        H3(1,lg0+k,lg0+k+j)=ms1;
                                        H3(2,lg0+k,lg0+k+j)=ms2;
                                        H3(1,lg0+k+j,lg0+k)=ms1;
                                        H3(2,lg0+k+j,lg0+k)=ms2;
                                    end
                                end
                            end    
                        end
                        ns = 0;
                        s1 = zeros(1,1);
                        s2 = zeros(1,1);
                        for j=1:lg0
                            for k=1:lg0
                                if (H3(2,lg0+j,lg0+k)>0)
                                    ns = ns + 1;
                                    s1(ns) = H3(1,lg0+j,lg0+k);
                                    s2(ns) = H3(2,lg0+j,lg0+k);
                                end
                            end
                        end
                        if (ns>0)
                            ms1 = median(s1);
                            ms2 = median(s2);
                            for j=1:lg0
                                for k=1:lg0
                                    if (H3(2,lg0+j,lg0+k)<=0)
                                        H3(1,lg0+j,lg0+k)=ms1;
                                        H3(2,lg0+j,lg0+k)=ms2;
                                    end
                                end
                            end
                        end
                        
                        
                        for j=1:lg0
                            ns = 0;
                            s1 = zeros(1,1);
                            s2 = zeros(1,1);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)>0)
                                    ns = ns + 1;
                                    s1(ns) = H3(1,k,j);
                                    s2(ns) = H3(2,k,j);
                                end
                            end
                            ms1 = median(s1);
                            ms2 = median(s2);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)==0)
                                    sm1 = (ms1/2)*(0.5*rand()+1.5);
                                    sm2 = (ms2/2)*(0.5*rand()+1.5);
                                    H3(1,k,j)=sm1;
                                    H3(1,j,k)=sm1;
                                    H3(2,k,j)=sm2;
                                    H3(2,j,k)=sm2;
                                end
                            end
                        end
                        for j=2*lg0+1:3*lg0
                            ns = 0;
                            s1 = zeros(1,1);
                            s2 = zeros(1,1);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)>0)
                                    ns = ns + 1;
                                    s1(ns) = H3(1,k,j);
                                    s2(ns) = H3(2,k,j);
                                end
                            end
                            ms1 = median(s1);
                            ms2 = median(s2);
                            for k=lg0+1:2*lg0
                                if (H3(2,k,j)==0)
                                    sm1 = (ms1/2)*(0.5*rand()+1.5);
                                    sm2 = (ms2/2)*(0.5*rand()+1.5);
                                    H3(1,k,j)=sm1;
                                    H3(1,j,k)=sm1;
                                    H3(2,k,j)=sm2;
                                    H3(2,j,k)=sm2;
                                end
                            end
                        end
                        
                        ngm = ngm + 1;
                        for j=1:3*lg0
                            for k=1:3*lg0
                                Hgmm(ngm,j,k) =  H3(1,j,k);
                                Hgm0m(ngm,j,k) =  H3(2,j,k);
                            end
                        end
                    end 
        end
    end
end

Hgt = zeros(2*lg1+lg0,2*lg1+lg0);
for i=1:2*lg1+lg0
    for j=1:2*lg1+lg0
        AAA = zeros(1,1);
        kk = 0;
        for k=1:ngp
            if(Hgp0m(k,i,j)>0)
                kk = kk + 1;
                AAA(kk) = Hgp0m(k,i,j);
            end
        end
        for k=1:ngm
            if(Hgm0m(k,2*lg1+lg0-i+1,2*lg1+lg0-j+1)>0)
                kk = kk + 1;
                AAA(kk) = Hgm0m(k,2*lg1+lg0-i+1,2*lg1+lg0-j+1);
            end
        end
        Hgt(i,j) = median(AAA);
    end
end

imagesc(log2(Hgt),[-1 1])


h = colorbar;
h.Label.String = 'log_2Hi-C(obs/exp)';
h.FontSize=15;
hold on
plot([lg1+1 lg1+1],[0 2*lg1+lg0],'k--','LineWidth',3)
plot([lg1+lg0 lg1+lg0],[0 2*lg1+lg0],'k--','LineWidth',3)
plot([0 2*lg1+lg0],[lg1+1 lg1+1],'k--','LineWidth',3)
plot([0 2*lg1+lg0],[lg1+lg0 lg1+lg0],'k--','LineWidth',3)
set (gca,'XTick',[lg1+1 lg1+lg0], 'XTickLabel' , {'TSS', 'TTS'});
set (gca,'YTick',[lg1+1 lg1+lg0], 'YTickLabel' , {'TSS', 'TTS'});
set(gca,'linewidth',3,'fontsize',20);


