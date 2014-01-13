#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;

require 'create_and_save_surface_and_charges.pm';
require 'read_settings.pm';

my $settings_file_name = $ARGV[0];
my %settings = &ReadSettings($settings_file_name);

my $document = DiscoveryScript::Open(
    {
        Path => $settings{"file_with_set"}      
    }
);

my $array_of_molecules = $document->Molecules;

for (my $ind = 0; $ind < $array_of_molecules->Count(); $ind++)
{
	$array_of_molecules->Item($ind)->SetInvisible();
}
# main cycle
for (my $ind = 0; $ind < $array_of_molecules->Count(); $ind++)
{
	my $molecule = $array_of_molecules->Item($ind);
	$settings{"curr_ind"} = $ind;
	&CreateAndSaveSurfaceAndCharges($document, $molecule, %settings);	
}