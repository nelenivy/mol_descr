set settings_file=E:\Disser\Tests\bzr\settings.txt
set scripts_folder=E:\Disser\mol_descr\DiscoveryStudioScripts\
set DS_perl_path="C:\Program Files (x86)\Accelrys\Discovery Studio 3.5\bin\perl.bat"
%DS_perl_path% "%scripts_folder%process_sdf.pl" %settings_file%
pause()