#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;
use POSIX;

sub GetMoleculeFullPath
{
#set_folder		- folder where set is 
#set_prefix		- prefix of the each molecule in set
#curr_ind 		- index of the molecule in the set
#set_size 		- set size
#digits_num_fixed 	- if true, number of digits is fixed and zeros are added for the missed digits
#digits_num 		- is used if digits_num is true
#separate_folder 	- is set true, if each file is situated in the separate folder
#extension		- file extension (mol, mol2, etc.)
	my (%settings, $curr_ind) = @_;
	my $set_folder = $settings{"set_folder"};
	my $set_prefix = $settings{"set_prefix"};
	my $extension = $settings{"extension"};
	my $digits_num_fixed = $settings{"digits_num_fixed"};
	my $digits_num = $settings{"digits_num"};
	my $separate_folder = $settings{"separate_folder"};
	
	if ((substr($set_folder, length($set_folder) - 1, 1) ne '\\') && (substr($set_folder, length($set_folder) - 1, 1) ne '/'))
	{
		$set_folder = $set_folder.'\\'
	}
	#create additional digits if need
	my $index_prefix = '';
	if ($digits_num_fixed)
	{
		my $zeros_needed = $digits_num - 1;
		
		if ($curr_ind > 0)
		{
			my $digits_in_index = floor(log($curr_ind) / log(10)) + 1;
			$zeros_needed = $digits_num - $digits_in_index;
		}
			
		for(my $i = 0; $i < $zeros_needed; $i++)
		{
		    $index_prefix = $index_prefix."0";
		}		
	}
	
	my $molecule_full_path = '';
	
	if ($separate_folder)
	{
		$molecule_full_path = $set_folder.$set_prefix.$index_prefix."$curr_ind".'\\'.$set_prefix.$index_prefix."$curr_ind".$extension;
	}
	else
	{
		$molecule_full_path = $set_folder.$set_prefix.$index_prefix."$curr_ind".$extension;
	}
	
	return $molecule_full_path;
}

1;