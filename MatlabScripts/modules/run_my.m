function run_my()


%OT1main('glik_1.txt','E:\Disser\Sets\Activity\glik110.sdf');
%OT2main('glik_2.txt','demo_distance.txt','E:\Disser\Sets\Activity\glik110.sdf');
%a= load('E:\Disser\mol_descr\MatlabScripts\modules\glik\MDmatrix\linear_fragments\level2\4\dbr\matrix.mat');
%b = cell2mat(struct2cell(a));
%save('E:\Disser\mol_descr\MatlabScripts\modules\glikmatrix1.txt','b','-ascii');

%OT1main('pirim_1.txt','E:\Disser\Sets\Activity\pirimidines205.sdf');
%OT2main('pirim_2.txt','demo_distance.txt','E:\Disser\Sets\Activity\pirimidines205.sdf');

%a1= load('E:\Disser\mol_descr\MatlabScripts\modules\pirim\MDmatrix\linear_fragments\level2\4\dbr\matrix.mat');
%b1 = cell2mat(struct2cell(a1));
%save('E:\Disser\mol_descr\MatlabScripts\modules\pirimmatrix.txt','b1','-ascii');

OT1main('bzr_1.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\bzr_3d.sd');
OT2main('bzr_2.txt','demo_distance.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\bzr_3d.sd');

a2= load('E:\Disser\mol_descr\MatlabScripts\modules\bzr\MDmatrix\linear_fragments\level2\4\dbr\matrix.mat');
b2 = cell2mat(struct2cell(a2));
save('E:\Disser\mol_descr\MatlabScripts\modules\bzrmatrix1.txt','b2','-ascii');

OT1main('er_lit_1.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\er_lit_3d.sd');
OT2main('er_lit_2.txt','demo_distance.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\er_lit_3d.sd');

a3= load('E:\Disser\mol_descr\MatlabScripts\modules\er_lit\MDmatrix\linear_fragments\level2\4\dbr\matrix.mat');
b3 = cell2mat(struct2cell(a3));
save('E:\Disser\mol_descr\MatlabScripts\modules\er_litmatrix1.txt','b3','-ascii');

OT1main('cox_1.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\cox2_3d.sd');
OT2main('cox_2.txt','demo_distance.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\cox2_3d.sd');

a4= load('E:\Disser\mol_descr\MatlabScripts\modules\cox\MDmatrix\linear_fragments\level2\4\dbr\matrix.mat');
b4 = cell2mat(struct2cell(a4));
save('E:\Disser\mol_descr\MatlabScripts\modules\coxmatrix1.txt','b4','-ascii');

OT1main('er_tox_1.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\er_tox_3d.sd');
OT2main('er_tox_2.txt','demo_distance.txt','E:\Disser\Sets\Activity\pharmocophore_kernel\er_tox_3d.sd');

a5= load('E:\Disser\mol_descr\MatlabScripts\modules\er_tox\MDmatrix\linear_fragments\level2\4\dbr\matrix.mat');
b5 = cell2mat(struct2cell(a5));
save('E:\Disser\mol_descr\MatlabScripts\modules\er_toxmatrix1.txt','b5','-ascii');

