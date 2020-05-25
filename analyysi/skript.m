load friction_data

%Laitetaan varied_force ja velocity data k‰tev‰sti taulukoihin.
%Taulukoissa vaihdeltu suure kasvaa indeksin mukana.

F = [[Favrforce9{:,1}]'...
    [Favrforce8{:,1}]'...
    [Favrforce7{:,1}]'...
    [Favrforce6{:,1}]'...
    [Favrforce5{:,1}]'...
    [Favrforce4{:,1}]'...
    [Favrforce3{:,1}]'...
    [Favrforce2{:,1}]'...
    [Favrforce1{:,1}]'...
    [Favrforce0{:,1}]'];

V = {[avrforce9{:,1}]'...
    [avrforce8{:,1}]'...
    [avrforce7{:,1}]'...
    [avrforce6{:,1}]'...
    [avrforce5{:,1}]'...
    [avrforce4{:,1}]'...
    [avrforce3{:,1}]'...
    [avrforce2{:,1}]'...
    [avrforce1{:,1}]'...
    [avrforce0{:,1}]'};

%Analysoidaan voimadata.


t = (1:size(F,1))*4*1e-3;
plots = plot(t,F);
legend([plots(1),plots(10)],{'smallest normal force /weigth of the slab','largest normal force /10x of that'})
xlabel('time /ps')
ylabel('frictional force /(eV/≈)')

%Etsit‰‰n kohta ja katkaistaan data.
F_search = sum(abs(F),2);
F_search = F_search(1e+4:1.2e+4);
[~, n] = min(F_search);
n = n + 1e+4;
F_cut = F(1:n,:);
t_cut = t(1:n);

figure(2)
plots2 = plot(t_cut,F_cut);
legend([plots2(1),plots2(10)],{'smallest normal force = weigth of the slab','largest normal force = 10x of that'})
xlabel('time /ps')
ylabel('frictional force /(eV/≈)')

%Keskiarvovoimat katkaistulle datalle.
F_avr_force = mean(F_cut);

figure(3)
plot(F_avr_force)
xlabel('normal force /weigth of the slab')
ylabel('average frictional force /(eV/≈)')



% Analysoidaan nopeusdata.

fig3 = figure(4);
hold on

V_cut = cell(10,1);
X = V_cut;
X_cut = V_cut;
V_avr_force = zeros(10,1);

for i = 1:10   
   f = V{i};
   
   %Etsit‰‰n kohta ja katkaistaan data.
   n = round(length(f)*(3/4));
   f_search = abs(f);
   f_search = f_search(n:end);
   [~, ind] = min(f_search);
   
   V_cut{i} = f(1:n+ind);
   
   %p‰‰llimm‰isten atomien paikka
   X{i} = (1:length(f))*i; %Hitaimman nopeuden timestepeiss‰ *(m/s).
   X{i} = X{i}*4/1e+4; %≈
   X_cut{i} = X{i}(1:n+ind);
   
   %Keskim‰‰r‰inen kitkavoima
   V_avr_force(i) = mean(V_cut{i});
   
   plots3(i) = plot(X_cut{i},V_cut{i});
   xlabel('x-coordinate of the top atoms /≈')
   ylabel('friction force /(eV/≈)')
end
legend([plots3(1),plots3(10)],{'smallest velocity = m/s','largest velocity = 10m/s'})

figure(5)
plot(V_avr_force)
xlabel('speed of the slab /(m/s)')
ylabel('average frictional force /(eV/≈)')



