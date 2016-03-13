###################################################################################
#: Title       : parseCHARMMpar 						  #
#: Author      : Salvatore Cosseddu - S.M.Cosseddu@warwick.ac.uk		  #
#: Institution : University of Warwick - Centre for Scientific Computing and	  #
#                                        School of Engineering			  #
#: Description : proc to parse CHARMM parameter file and load it values in the    #
#                VMD environment                                                  #
#: Date Created: Tue Dec 13 17:49:58 GMT 2011					  #
#: Last Edit   : Tue Dec 13 17:49:58 GMT 2011    				  #
#: Version     : 0.01								  #
# 										  #
#										  #
#  COPYRIGHT					       				  #
#  Copyright © 2010 Salvatore Cosseddu		       				  #
#  Centre for Scientific Computing, University of Warwick.		       	  #
#  License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>. #
#  This is free software: you are free to change and redistribute it.         	  #
#  There is NO WARRANTY, to the extent permitted by law.          		  #
#										  #
###################################################################################

namespace eval ::parseCHARMMpar:: {
    variable parameterType
    variable lj_parameter
    array set parameterType {
	bonds 0
	angles 0
	dihedrals 0
	improper 0
	cmap 0
	lj 0 
	nbfix 0 
    } 
    
    proc zeroingArray {workarray} {
	foreach element [array name $workarray] {
	    set $workarray($element) 0
	}
    }
}


########################################
##       parsing file                 ##
########################################
proc ::parseCHARMMpar::parsefile { parafile } {

    variable parameterType
    variable lj_parameter

    # open para file
    if { [ catch {open $parafile r} infile ] } { error "error opening $parafile" }
    
    puts "Reading from $parafile"
    # reading file
    foreach line [split [read $infile] \n] {
	
	set line [string trim $line]
	
	# Skip comments
	if {[string index $line 0]=="!"} { continue }
	if {[string index $line 0]=="*"} { continue }

	set keyword [string toupper [lindex $line 0]]
	# checking the keyword present in the parameter file.
	# 
	switch -exact -- $keyword {
	    "BONDS"       { 
		zeroingArray parameterType
		set parameterType(bonds) 1
	    }
	    "ANGLES"	  { 
		zeroingArray parameterType
		set parameterType(angles) 1
	    }
	    "DIHEDRALS"	  { 
		zeroingArray parameterType
		set parameterType(dihedrals) 1
	    }
	    "IMPROPER"	  { 
		zeroingArray parameterType
		set parameterType(improper) 1
	    }
	    "CMAP"	  { 
		zeroingArray parameterType
		set parameterType(cmap) 1
	    }
	    "NONBONDED"   { 
		zeroingArray parameterType
		set parameterType(lj) 1
	    }
	    "NBFIX"	  { 
		zeroingArray parameterType
		set parameterType(nbfix) 1
	    }
	}
	
	# if NONBONDED keyword is found read eps Rmin/2
	if {$parameterType(lj)} {
	    # parsing lj data
	    if {[ regexp {^[a-zA-Z0-9]+\s+([+-]*[0-9]+\.[0-9]+\s*){3}} $line ]} {
		lassign $line atom_type ignored epsilon	Rmin   
		puts "Read LJ par for $atom_type"
		set lj_parameter($atom_type) [list $epsilon $Rmin]
	    }
	}
    }

}

proc ::parseCHARMMpar::VdWassign {molID parafile} {
    variable parameterType
    variable lj_parameter

    puts "parsing $parafile"
    ::parseCHARMMpar::parsefile $parafile
    puts "done..."
    set all [atomselect $molID all]
    
    set typelist [lsort -unique [$all get type]]
    puts "the following atom type are in the molecule $molID
  $typelist"

    foreach type $typelist {
	set tmp [atomselect $molID "type $type"]
	set Rmin [lindex $lj_parameter($type) 1]
	$tmp set radius $Rmin
	$tmp delete
    }
    $all delete
    puts "done"
}