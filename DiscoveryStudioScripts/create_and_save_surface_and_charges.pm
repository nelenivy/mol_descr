#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;

sub CreateAndSaveSurfaceChargesAndWDVRadii
{
#document - document which contains molecule
	my ($document, $molecule, %settings) = @_;
	my $curr_ind = $settings{"curr_ind"};
	$molecule->SetVisible();
	my $output_folder = $settings{"output_folder"};
	my $prefix = $settings{"output_prefix"};
	#add slash if need
	if ((substr($output_folder, length($output_folder) - 1, 1) ne '\\') && (substr($output_folder, length($output_folder) - 1, 1) ne '/'))
	{
		$output_folder = $output_folder.'\\'
	}
	#surface
	my $array_of_atoms = $molecule->Atoms; 
	my $surface = $document->CreateSolidSurface( $array_of_atoms, Mdm::surfaceStyleSolvent, False,
	    Mdm::surfaceColorByElectrostaticPotential, 1.4 );	
	$document->Save($output_folder.$prefix."_$curr_ind".'.wrl' ,'wrl' ); 
	#charges
	$document->CalculateCharges();
	open(my $charges_file, '>', $output_folder.$prefix."_$curr_ind".'.ch') ;

	foreach my $atom (@$array_of_atoms)
	{
	    print $charges_file $atom->XYZ->X.' '.$atom->XYZ->Y.' '.$atom->XYZ->Z.' '.$atom->Charge."\n";
	}
	
	close $charges_file;
	#wdv radii
	open(my $radii_file, '>', $output_folder.$prefix."_$curr_ind".'.wdv') ;

	foreach my $atom (@$array_of_atoms)
	{
	    print $radii_file $atom->XYZ->X.' '.$atom->XYZ->Y.' '.$atom->XYZ->Z.' '.$atom->VdwRadius."\n";
	}
	
	close $radii_file;
	
	$surface->SetInvisible();
	$molecule->SetInvisible();
	$surface = undef;
	$molecule = undef;
	$array_of_atoms = undef;
	$document = undef;
	return 1;
}

1;