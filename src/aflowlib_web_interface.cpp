// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// fixed for XZ - SC 2018-2019
// fixed for new AUID language (SC 2019)
// fixed for tree search on the AUID directories (SC 2019) super-speed

#ifndef _AFLOWLIB_WEB_INTERFACE_CPP_
#define _AFLOWLIB_WEB_INTERFACE_CPP_
#include "aflow.h"
#include "aflowlib_webapp_entry.cpp"  //CO20170622 - BH JMOL stuff
#include "aflowlib_webapp_bands.cpp"  //CO20180305 - GG bands stuff

const string _DEVIL_PROTOTYPES_STRING_ = "64,65,549,550,f8269,f9083,f8819";

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ***************************************************************************
namespace aflowlib {
  //  class _aflowlib_entry 

  _aflowlib_entry::_aflowlib_entry() {  // constructor PUBLIC
    entry.clear();ventry.clear();
    auid.clear();
    vauid.clear();vauid.clear();
    aurl.clear();vaurl.clear();
    keywords.clear();vkeywords.clear();
    aflowlib_date.clear();vaflowlib_date.clear(); //CO20200624 - adding LOCK date
    aflowlib_version.clear();
    aflowlib_entries.clear();vaflowlib_entries.clear();
    aflowlib_entries_number=0;
    aflow_version.clear();
    catalog.clear();
    data_api="aapi1.2"; // new version of the API
    data_source="aflowlib";
    data_language="";
    error_status.clear();
    author.clear();vauthor.clear();
    calculation_cores=1;
    calculation_memory=AUROSTD_NAN;
    calculation_time=AUROSTD_NAN;
    corresponding.clear();vcorresponding.clear();
    loop.clear();vloop.clear();
    node_CPU_Cores=AUROSTD_NAN;node_CPU_MHz=AUROSTD_NAN;node_CPU_Model.clear();node_RAM_GB=AUROSTD_NAN;
    Bravais_lattice_orig.clear();Bravais_lattice_relax.clear();
    code.clear();
    composition.clear();vcomposition.clear();
    compound.clear();
    density=AUROSTD_NAN;
    density_orig=AUROSTD_NAN; //DX20190124 - add original crystal info
    dft_type.clear();vdft_type.clear();
    eentropy_cell=AUROSTD_NAN;eentropy_atom=AUROSTD_NAN;
    Egap=AUROSTD_NAN;Egap_fit=AUROSTD_NAN;
    energy_cell=AUROSTD_NAN;energy_atom=AUROSTD_NAN;energy_atom_relax1=AUROSTD_NAN;
    energy_cutoff=AUROSTD_NAN;
    delta_electronic_energy_convergence=AUROSTD_NAN;
    delta_electronic_energy_threshold=AUROSTD_NAN;
    nkpoints=0;
    nkpoints_irreducible=0;
    kppra=0;
    kpoints.clear();
    kpoints_nnn_relax.clear();
    kpoints_nnn_static.clear();
    kpoints_pairs.clear();
    kpoints_bands_path_grid=0;
    enthalpy_cell=AUROSTD_NAN;enthalpy_atom=AUROSTD_NAN;
    enthalpy_formation_cell=AUROSTD_NAN;enthalpy_formation_atom=AUROSTD_NAN;
    entropic_temperature=AUROSTD_NAN;
    files.clear();vfiles.clear();
    files_LIB.clear();vfiles_LIB.clear();
    files_RAW.clear();vfiles_RAW.clear();
    files_WEB.clear();vfiles_WEB.clear();
    forces.clear();vforces.clear();
    Egap_type.clear();
    geometry.clear();vgeometry.clear();
    geometry_orig.clear();vgeometry_orig.clear(); //DX20190124 - add original crystal info
    lattice_system_orig.clear();lattice_variation_orig.clear();lattice_system_relax.clear();lattice_variation_relax.clear();
    ldau_TLUJ.clear();
    vLDAU.resize(4);  //ME20190129
    natoms=AUROSTD_NAN;
    natoms_orig=AUROSTD_NAN; //DX20190124 - add original crystal info
    nbondxx.clear();vnbondxx.clear();
    nspecies=AUROSTD_NAN;
    Pearson_symbol_orig.clear();Pearson_symbol_relax.clear();
    positions_cartesian.clear();vpositions_cartesian.clear();
    positions_fractional.clear();vpositions_fractional.clear();
    pressure=AUROSTD_NAN;
    stress_tensor.clear();vstress_tensor.clear();
    pressure_residual=AUROSTD_NAN;
    Pulay_stress=AUROSTD_NAN;
    prototype.clear();
    PV_cell=AUROSTD_NAN;PV_atom=AUROSTD_NAN;
    scintillation_attenuation_length=AUROSTD_NAN;
    sg.clear();sg2.clear();vsg.clear();vsg2.clear();  //CO20171202
    spacegroup_orig.clear();spacegroup_relax.clear();
    species.clear();vspecies.clear();
    species_pp.clear();vspecies_pp.clear();
    species_pp_version.clear();vspecies_pp_version.clear();
    species_pp_ZVAL.clear();vspecies_pp_ZVAL.clear();
    species_pp_AUID.clear();vspecies_pp_AUID.clear();
    METAGGA.clear();
    spin_cell=AUROSTD_NAN;spin_atom=AUROSTD_NAN;
    spinD.clear();vspinD.clear();
    spinD_magmom_orig.clear();vspinD_magmom_orig.clear();
    spinF=AUROSTD_NAN;
    sponsor.clear();vsponsor.clear();
    stoichiometry.clear();vstoichiometry.clear();
    valence_cell_std=AUROSTD_NAN;valence_cell_iupac=AUROSTD_NAN;
    volume_cell=AUROSTD_NAN;volume_atom=AUROSTD_NAN;
    volume_cell_orig=AUROSTD_NAN;volume_atom_orig=AUROSTD_NAN; //DX20190124 - add original crystal info
    //DX20190124 - added original symmetry info - START
    // SYMMETRY
    crystal_family_orig="";
    crystal_system_orig="";
    crystal_class_orig="";
    point_group_Hermann_Mauguin_orig="";
    point_group_Schoenflies_orig="";
    point_group_orbifold_orig="";
    point_group_type_orig="";
    point_group_order_orig=AUROSTD_NAN;
    point_group_structure_orig="";
    Bravais_lattice_lattice_type_orig="";
    Bravais_lattice_lattice_variation_type_orig="";
    Bravais_lattice_lattice_system_orig="";
    Bravais_superlattice_lattice_type_orig="";
    Bravais_superlattice_lattice_variation_type_orig="";
    Bravais_superlattice_lattice_system_orig="";
    Pearson_symbol_superlattice_orig="";
    reciprocal_lattice_type_orig="";
    reciprocal_lattice_variation_type_orig="";
    reciprocal_geometry_orig.clear();vreciprocal_geometry_orig.clear();
    reciprocal_volume_cell_orig=AUROSTD_NAN;
    Wyckoff_letters_orig="";
    Wyckoff_multiplicities_orig="";
    Wyckoff_site_symmetries_orig="";
    //DX20190124 - added original symmetry info - END
    //DX20180823 - added more symmetry info - START
    // SYMMETRY
    crystal_family="";
    crystal_system="";
    crystal_class="";
    point_group_Hermann_Mauguin="";
    point_group_Schoenflies="";
    point_group_orbifold="";
    point_group_type="";
    point_group_order=AUROSTD_NAN;
    point_group_structure="";
    Bravais_lattice_lattice_type="";
    Bravais_lattice_lattice_variation_type="";
    Bravais_lattice_lattice_system="";
    Bravais_superlattice_lattice_type="";
    Bravais_superlattice_lattice_variation_type="";
    Bravais_superlattice_lattice_system="";
    Pearson_symbol_superlattice="";
    reciprocal_lattice_type="";
    reciprocal_lattice_variation_type="";
    reciprocal_geometry.clear();vreciprocal_geometry.clear();
    reciprocal_volume_cell=AUROSTD_NAN;
    Wyckoff_letters="";
    Wyckoff_multiplicities="";
    Wyckoff_site_symmetries="";
    //DX20180823 - added more symmetry info - END
    //DX20190209 - added anrl info - START
    anrl_label_orig="";
    anrl_parameter_list_orig="";
    anrl_parameter_values_orig="";
    anrl_label_relax="";
    anrl_parameter_list_relax="";
    anrl_parameter_values_relax="";
    //DX20190209 - added anrl info - END
    // AGL/AEL
    agl_thermal_conductivity_300K=AUROSTD_NAN;
    agl_debye=AUROSTD_NAN;
    agl_acoustic_debye=AUROSTD_NAN;
    agl_gruneisen=AUROSTD_NAN;
    agl_heat_capacity_Cv_300K=AUROSTD_NAN;
    agl_heat_capacity_Cp_300K=AUROSTD_NAN;
    agl_thermal_expansion_300K=AUROSTD_NAN;
    agl_bulk_modulus_static_300K=AUROSTD_NAN;
    agl_bulk_modulus_isothermal_300K=AUROSTD_NAN;
    agl_poisson_ratio_source=""; //CT20181212
    agl_vibrational_free_energy_300K_cell=AUROSTD_NAN; //CT20181212
    agl_vibrational_free_energy_300K_atom=AUROSTD_NAN; //CT20181212
    agl_vibrational_entropy_300K_cell=AUROSTD_NAN; //CT20181212
    agl_vibrational_entropy_300K_atom=AUROSTD_NAN; //CT20181212
    ael_poisson_ratio=AUROSTD_NAN;
    ael_bulk_modulus_voigt=AUROSTD_NAN;
    ael_bulk_modulus_reuss=AUROSTD_NAN;
    ael_shear_modulus_voigt=AUROSTD_NAN;
    ael_shear_modulus_reuss=AUROSTD_NAN;
    ael_bulk_modulus_vrh=AUROSTD_NAN;
    ael_shear_modulus_vrh=AUROSTD_NAN;
    ael_elastic_anisotropy=AUROSTD_NAN; //CO20181129
    ael_youngs_modulus_vrh=AUROSTD_NAN; //CT20181212
    ael_speed_sound_transverse=AUROSTD_NAN; //CT20181212
    ael_speed_sound_longitudinal=AUROSTD_NAN; //CT20181212
    ael_speed_sound_average=AUROSTD_NAN; //CT20181212
    ael_pughs_modulus_ratio=AUROSTD_NAN; //CT20181212
    ael_debye_temperature=AUROSTD_NAN; //CT20181212
    ael_applied_pressure=AUROSTD_NAN; //CT20181212
    ael_average_external_pressure=AUROSTD_NAN; //CT20181212
    ael_stiffness_tensor.clear();  //ME20191105
    ael_compliance_tensor.clear();  //ME20191105
    // BADER
    bader_net_charges.clear();vbader_net_charges.clear();
    bader_atomic_volumes.clear();vbader_atomic_volumes.clear();
    // legacy
    server.clear();vserver.clear();vserverdir.clear();
    icsd.clear();
    stoich.clear();vstoich.clear();
    // apennsy
    structure_name.clear();  // apennsy
    structure_description.clear();  // apennsy
    distance_gnd=AUROSTD_NAN;  // apennsy
    distance_tie=AUROSTD_NAN;  // apennsy
    pureA=FALSE;pureB=FALSE;  // apennsy
    fcc=FALSE; bcc=FALSE;hcp=FALSE;  // apennsy
    stoich_a=AUROSTD_NAN;stoich_b=AUROSTD_NAN;  // apennsy
    bond_aa=AUROSTD_NAN;bond_ab=AUROSTD_NAN;bond_bb=AUROSTD_NAN;  // apennsy
    vNsgroup.clear();  // apennsy
    vsgroup.clear();  // apennsy
    vstr.clear();  // apennsy
  }
  _aflowlib_entry::~_aflowlib_entry() { // destructor PUBLIC
    free();
  }

  void _aflowlib_entry::copy(const _aflowlib_entry& b) { // copy PRIVATE
    entry=b.entry;ventry.clear();for(uint i=0;i<b.ventry.size();i++) ventry.push_back(b.ventry.at(i));
    auid=b.auid;
    vauid=b.vauid;vauid.clear();for(uint i=0;i<b.vauid.size();i++) vauid.push_back(b.vauid.at(i));
    aurl=b.aurl;vaurl.clear();for(uint i=0;i<b.vaurl.size();i++) vaurl.push_back(b.vaurl.at(i));
    vkeywords.clear();for(uint i=0;i<b.vkeywords.size();i++) vkeywords.push_back(b.vkeywords.at(i));
    aflowlib_date=b.aflowlib_date;for(uint i=0;i<b.vaflowlib_date.size();i++) vaflowlib_date.push_back(b.vaflowlib_date.at(i)); //CO20200624 - adding LOCK date
    aflowlib_version=b.aflowlib_version;
    aflowlib_entries=b.aflowlib_entries;
    vaflowlib_entries.clear();for(uint i=0;i<b.vaflowlib_entries.size();i++) vaflowlib_entries.push_back(b.vaflowlib_entries.at(i));
    aflowlib_entries_number=b.aflowlib_entries_number;
    aflow_version=b.aflow_version;
    catalog=b.catalog;
    data_api=b.data_api;
    data_source=b.data_source;
    data_language=b.data_language;
    error_status=b.error_status;
    author=b.author;vauthor.clear();for(uint i=0;i<b.vauthor.size();i++) vauthor.push_back(b.vauthor.at(i));
    calculation_cores=b.calculation_cores;calculation_memory=b.calculation_memory;calculation_time=b.calculation_time;
    corresponding=b.corresponding;vcorresponding.clear();for(uint i=0;i<b.vcorresponding.size();i++) vcorresponding.push_back(b.vcorresponding.at(i));
    loop=b.loop;vloop.clear();for(uint i=0;i<b.vloop.size();i++) vloop.push_back(b.vloop.at(i));
    node_CPU_Cores=b.node_CPU_Cores;node_CPU_MHz=b.node_CPU_MHz;node_CPU_Model=b.node_CPU_Model;node_RAM_GB=b.node_RAM_GB;
    Bravais_lattice_orig=b.Bravais_lattice_orig;Bravais_lattice_relax=b.Bravais_lattice_relax;
    code=b.code;
    composition=b.composition;vcomposition.clear();for(uint i=0;i<b.vcomposition.size();i++) vcomposition.push_back(b.vcomposition.at(i));
    compound=b.compound;
    density=b.density;
    density_orig=b.density_orig; //DX20190124 - add original crystal info
    dft_type=b.dft_type;vdft_type.clear();for(uint i=0;i<b.vdft_type.size();i++) vdft_type.push_back(b.vdft_type.at(i));
    eentropy_cell=b.eentropy_cell;eentropy_atom=b.eentropy_atom;
    Egap=b.Egap;Egap_fit=b.Egap_fit;
    energy_cell=b.energy_cell;energy_atom=b.energy_atom;energy_atom_relax1=b.energy_atom_relax1;
    energy_cutoff=b.energy_cutoff;
    delta_electronic_energy_convergence=b.delta_electronic_energy_convergence;
    delta_electronic_energy_threshold=b.delta_electronic_energy_threshold;
    nkpoints=b.nkpoints;
    nkpoints_irreducible=b.nkpoints_irreducible;
    kppra=b.kppra;
    kpoints=b.kpoints;
    kpoints_nnn_relax=b.kpoints_nnn_relax;
    kpoints_nnn_static=b.kpoints_nnn_static;
    kpoints_pairs=b.kpoints_pairs;
    kpoints_bands_path_grid=b.kpoints_bands_path_grid;
    enthalpy_cell=b.enthalpy_cell;enthalpy_atom=b.enthalpy_atom;
    enthalpy_formation_cell=b.enthalpy_formation_cell;enthalpy_formation_atom=b.enthalpy_formation_atom;
    entropic_temperature=b.entropic_temperature;
    files=b.files;vfiles.clear();for(uint i=0;i<b.vfiles.size();i++) vfiles.push_back(b.vfiles.at(i));
    files_LIB=b.files_LIB;vfiles_LIB.clear();for(uint i=0;i<b.vfiles_LIB.size();i++) vfiles_LIB.push_back(b.vfiles_LIB.at(i));
    files_RAW=b.files_RAW;vfiles_RAW.clear();for(uint i=0;i<b.vfiles_RAW.size();i++) vfiles_RAW.push_back(b.vfiles_RAW.at(i));
    files_WEB=b.files_WEB;vfiles_WEB.clear();for(uint i=0;i<b.vfiles_WEB.size();i++) vfiles_WEB.push_back(b.vfiles_WEB.at(i));
    forces=b.forces;vforces.clear();for(uint i=0;i<b.vforces.size();i++) vforces.push_back(b.vforces.at(i));
    Egap_type=b.Egap_type;
    geometry=b.geometry;vgeometry.clear();for(uint i=0;i<b.vgeometry.size();i++) vgeometry.push_back(b.vgeometry.at(i));
    geometry_orig=b.geometry_orig;vgeometry_orig.clear();for(uint i=0;i<b.vgeometry_orig.size();i++) vgeometry_orig.push_back(b.vgeometry_orig.at(i)); //DX20190124 - add original crystal info
    lattice_system_orig=b.lattice_system_orig;lattice_variation_orig=b.lattice_variation_orig;
    lattice_system_relax=b.lattice_system_relax;lattice_variation_relax=b.lattice_variation_relax;
    ldau_TLUJ=b.ldau_TLUJ;
    vLDAU=b.vLDAU;  //ME20190129
    natoms=b.natoms;
    natoms_orig=b.natoms_orig; //DX20190124 - add original crystal info
    nbondxx=b.nbondxx;vnbondxx.clear();for(uint i=0;i<b.vnbondxx.size();i++) vnbondxx.push_back(b.vnbondxx.at(i));
    nspecies=b.nspecies;
    Pearson_symbol_orig=b.Pearson_symbol_orig;Pearson_symbol_relax=b.Pearson_symbol_relax;
    positions_cartesian=b.positions_cartesian;vpositions_cartesian.clear();for(uint i=0;i<b.vpositions_cartesian.size();i++) vpositions_cartesian.push_back(b.vpositions_cartesian.at(i));
    positions_fractional=b.positions_fractional;vpositions_fractional.clear();for(uint i=0;i<b.vpositions_fractional.size();i++) vpositions_fractional.push_back(b.vpositions_fractional.at(i));
    pressure=b.pressure;
    stress_tensor=b.stress_tensor;vstress_tensor.clear();for(uint i=0;i<b.vstress_tensor.size();i++) vstress_tensor.push_back(b.vstress_tensor.at(i));
    pressure_residual=b.pressure_residual;
    Pulay_stress=b.Pulay_stress;
    prototype=b.prototype;
    PV_cell=b.PV_cell;PV_atom=b.PV_atom;
    scintillation_attenuation_length=b.scintillation_attenuation_length;
    sg=b.sg;sg2=b.sg2;vsg.clear();for(uint i=0;i<b.vsg.size();i++){vsg.push_back(b.vsg[i]);} vsg2.clear();for(uint i=0;i<b.vsg2.size();i++){vsg2.push_back(b.vsg2[i]);}  //CO20171202
    spacegroup_orig=b.spacegroup_orig;spacegroup_relax=b.spacegroup_relax;
    species=b.species;vspecies.clear();for(uint i=0;i<b.vspecies.size();i++) vspecies.push_back(b.vspecies.at(i));
    species_pp=b.species_pp;vspecies_pp.clear();for(uint i=0;i<b.vspecies_pp.size();i++) vspecies_pp.push_back(b.vspecies_pp.at(i));
    species_pp_version=b.species_pp_version;vspecies_pp_version.clear();for(uint i=0;i<b.vspecies_pp_version.size();i++) vspecies_pp_version.push_back(b.vspecies_pp_version.at(i));
    species_pp_ZVAL=b.species_pp_ZVAL;vspecies_pp_ZVAL.clear();for(uint i=0;i<b.vspecies_pp_ZVAL.size();i++) vspecies_pp_ZVAL.push_back(b.vspecies_pp_ZVAL.at(i));
    species_pp_AUID=b.species_pp_AUID;vspecies_pp_AUID.clear();for(uint i=0;i<b.vspecies_pp_AUID.size();i++) vspecies_pp_AUID.push_back(b.vspecies_pp_AUID.at(i));
    METAGGA=b.METAGGA;
    spin_cell=b.spin_cell;spin_atom=b.spin_atom;
    spinD=b.spinD;vspinD.clear();for(uint i=0;i<b.vspinD.size();i++) vspinD.push_back(b.vspinD.at(i));
    spinD_magmom_orig=b.spinD_magmom_orig;vspinD_magmom_orig.clear();for(uint i=0;i<b.vspinD_magmom_orig.size();i++) vspinD_magmom_orig.push_back(b.vspinD_magmom_orig.at(i));
    spinF=b.spinF;
    sponsor=b.sponsor;vsponsor.clear();for(uint i=0;i<b.vsponsor.size();i++) vsponsor.push_back(b.vsponsor.at(i));
    stoichiometry=b.stoichiometry;vstoichiometry.clear();for(uint i=0;i<b.vstoichiometry.size();i++) vstoichiometry.push_back(b.vstoichiometry.at(i));
    valence_cell_std=b.valence_cell_std;valence_cell_iupac=b.valence_cell_iupac;
    volume_cell=b.volume_cell;volume_atom=b.volume_atom;
    volume_cell_orig=b.volume_cell_orig;volume_atom_orig=b.volume_atom_orig; //DX20190124 - add original crystal info
    //DX20190124 - added original symmetry info - START
    // SYMMETRY
    crystal_family_orig=b.crystal_family_orig;
    crystal_system_orig=b.crystal_system_orig;
    crystal_class_orig=b.crystal_class_orig;
    point_group_Hermann_Mauguin_orig=b.point_group_Hermann_Mauguin_orig;
    point_group_Schoenflies_orig=b.point_group_Schoenflies_orig;
    point_group_orbifold_orig=b.point_group_orbifold_orig;
    point_group_type_orig=b.point_group_type_orig;
    point_group_order_orig=b.point_group_order_orig;
    point_group_structure_orig=b.point_group_structure_orig;
    Bravais_lattice_lattice_type_orig=b.Bravais_lattice_lattice_type_orig;
    Bravais_lattice_lattice_variation_type_orig=b.Bravais_lattice_lattice_variation_type_orig;
    Bravais_lattice_lattice_system_orig=b.Bravais_lattice_lattice_system_orig;
    Bravais_superlattice_lattice_type_orig=b.Bravais_superlattice_lattice_type_orig;
    Bravais_superlattice_lattice_variation_type_orig=b.Bravais_superlattice_lattice_variation_type_orig;
    Bravais_superlattice_lattice_system_orig=b.Bravais_superlattice_lattice_system_orig;
    Pearson_symbol_superlattice_orig=b.Pearson_symbol_superlattice_orig;
    reciprocal_geometry_orig=b.reciprocal_geometry_orig;vreciprocal_geometry_orig.clear();for(uint i=0;i<b.vreciprocal_geometry_orig.size();i++) vreciprocal_geometry_orig.push_back(b.vreciprocal_geometry_orig.at(i));
    reciprocal_volume_cell_orig=b.reciprocal_volume_cell_orig;
    reciprocal_lattice_type_orig=b.reciprocal_lattice_type_orig;
    reciprocal_lattice_variation_type_orig=b.reciprocal_lattice_variation_type_orig;
    Wyckoff_letters_orig=b.Wyckoff_letters_orig;
    Wyckoff_multiplicities_orig=b.Wyckoff_multiplicities_orig;
    Wyckoff_site_symmetries_orig=b.Wyckoff_site_symmetries_orig;
    //DX20190124 - added original symmetry info - END
    //DX20180823 - added more symmetry info - START
    // SYMMETRY
    crystal_family=b.crystal_family;
    crystal_system=b.crystal_system;
    crystal_class=b.crystal_class;
    point_group_Hermann_Mauguin=b.point_group_Hermann_Mauguin;
    point_group_Schoenflies=b.point_group_Schoenflies;
    point_group_orbifold=b.point_group_orbifold;
    point_group_type=b.point_group_type;
    point_group_order=b.point_group_order;
    point_group_structure=b.point_group_structure;
    Bravais_lattice_lattice_type=b.Bravais_lattice_lattice_type;
    Bravais_lattice_lattice_variation_type=b.Bravais_lattice_lattice_variation_type;
    Bravais_lattice_lattice_system=b.Bravais_lattice_lattice_system;
    Bravais_superlattice_lattice_type=b.Bravais_superlattice_lattice_type;
    Bravais_superlattice_lattice_variation_type=b.Bravais_superlattice_lattice_variation_type;
    Bravais_superlattice_lattice_system=b.Bravais_superlattice_lattice_system;
    Pearson_symbol_superlattice=b.Pearson_symbol_superlattice;
    reciprocal_geometry=b.reciprocal_geometry;vreciprocal_geometry.clear();for(uint i=0;i<b.vreciprocal_geometry.size();i++) vreciprocal_geometry.push_back(b.vreciprocal_geometry.at(i));
    reciprocal_volume_cell=b.reciprocal_volume_cell; //DX20190124 - fix typo, add reciprocal
    reciprocal_lattice_type=b.reciprocal_lattice_type;
    reciprocal_lattice_variation_type=b.reciprocal_lattice_variation_type;
    Wyckoff_letters=b.Wyckoff_letters;
    Wyckoff_multiplicities=b.Wyckoff_multiplicities;
    Wyckoff_site_symmetries=b.Wyckoff_site_symmetries;
    //DX20180823 - added more symmetry info - END
    //DX20190209 - added anrl info - START
    anrl_label_orig=b.anrl_label_orig;
    anrl_parameter_list_orig=b.anrl_parameter_list_orig;
    anrl_parameter_values_orig=b.anrl_parameter_values_orig;
    anrl_label_relax=b.anrl_label_relax;
    anrl_parameter_list_relax=b.anrl_parameter_list_relax;
    anrl_parameter_values_relax=b.anrl_parameter_values_relax;
    //DX20190209 - added anrl info - END
    // AGL/AEL
    agl_thermal_conductivity_300K=b.agl_thermal_conductivity_300K;
    agl_debye=b.agl_debye;
    agl_acoustic_debye=b.agl_acoustic_debye;
    agl_gruneisen=b.agl_gruneisen;
    agl_heat_capacity_Cv_300K=b.agl_heat_capacity_Cv_300K;
    agl_heat_capacity_Cp_300K=b.agl_heat_capacity_Cp_300K;
    agl_thermal_expansion_300K=b.agl_thermal_expansion_300K;
    agl_bulk_modulus_static_300K=b.agl_bulk_modulus_static_300K;
    agl_bulk_modulus_isothermal_300K=b.agl_bulk_modulus_isothermal_300K;
    agl_poisson_ratio_source=b.agl_poisson_ratio_source; //CT20181212
    agl_vibrational_free_energy_300K_cell=b.agl_vibrational_free_energy_300K_cell; //CT20181212
    agl_vibrational_free_energy_300K_atom=b.agl_vibrational_free_energy_300K_atom; //CT20181212
    agl_vibrational_entropy_300K_cell=b.agl_vibrational_entropy_300K_cell; //CT20181212
    agl_vibrational_entropy_300K_atom=b.agl_vibrational_entropy_300K_atom; //CT20181212
    ael_poisson_ratio=b.ael_poisson_ratio;
    ael_bulk_modulus_voigt=b.ael_bulk_modulus_voigt;
    ael_bulk_modulus_reuss=b.ael_bulk_modulus_reuss;
    ael_shear_modulus_voigt=b.ael_shear_modulus_voigt;
    ael_shear_modulus_reuss=b.ael_shear_modulus_reuss;
    ael_bulk_modulus_vrh=b.ael_bulk_modulus_vrh;
    ael_shear_modulus_vrh=b.ael_shear_modulus_vrh;
    ael_elastic_anisotropy=b.ael_elastic_anisotropy; //CO20181129
    ael_youngs_modulus_vrh=b.ael_youngs_modulus_vrh; //CT20181212
    ael_speed_sound_transverse=b.ael_speed_sound_transverse; //CT20181212
    ael_speed_sound_longitudinal=b.ael_speed_sound_longitudinal; //CT20181212
    ael_speed_sound_average=b.ael_speed_sound_average; //CT20181212
    ael_pughs_modulus_ratio=b.ael_pughs_modulus_ratio; //CT20181212
    ael_debye_temperature=b.ael_debye_temperature; //CT20181212
    ael_applied_pressure=b.ael_applied_pressure; //CT20181212
    ael_average_external_pressure=b.ael_average_external_pressure; //CT20181212
    ael_stiffness_tensor = b.ael_stiffness_tensor;  //ME20191105
    ael_compliance_tensor = b.ael_compliance_tensor;  //ME20191105
    // BADER
    bader_net_charges=b.bader_net_charges;vbader_net_charges.clear();for(uint i=0;i<b.vbader_net_charges.size();i++) vbader_net_charges.push_back(b.vbader_net_charges.at(i));
    bader_atomic_volumes=b.bader_atomic_volumes;vbader_atomic_volumes.clear();for(uint i=0;i<b.vbader_atomic_volumes.size();i++) vbader_atomic_volumes.push_back(b.vbader_atomic_volumes.at(i));
    // legacy
    server=b.server;
    vserver.clear();for(uint i=0;i<b.vserver.size();i++) vserver.push_back(b.vserver.at(i));
    vserverdir.clear();for(uint i=0;i<b.vserverdir.size();i++) vserverdir.push_back(b.vserverdir.at(i));
    icsd=b.icsd;
    stoich=b.stoich;vstoich.clear();for(uint i=0;i<b.vstoich.size();i++) vstoich.push_back(b.vstoich.at(i));
    // apennsy
    structure_name=b.structure_name;  // apennsy
    structure_description=b.structure_description;  // apennsy
    distance_gnd=b.distance_gnd;  // apennsy
    distance_tie=b.distance_tie;  // apennsy
    pureA=b.pureA;pureB=b.pureB;  // apennsy
    fcc=b.fcc;bcc=b.bcc;hcp=b.hcp;  // apennsy
    stoich_a=b.stoich_a;stoich_b=b.stoich_b;  // apennsy
    bond_aa=b.bond_aa;bond_ab=b.bond_ab;bond_bb=b.bond_bb;  // apennsy
    vNsgroup.clear();for(uint i=0;i<b.vNsgroup.size();i++) vNsgroup.push_back(b.vNsgroup.at(i));  // apennsy
    vsgroup.clear();for(uint i=0;i<b.vsgroup.size();i++) vsgroup.push_back(b.vsgroup.at(i));  // apennsy
    vstr.clear();for(uint i=0;i<b.vstr.size();i++) vstr.push_back(b.vstr.at(i));  // apennsy
  }


  const _aflowlib_entry& _aflowlib_entry::operator=(const _aflowlib_entry& b) {  // operator= PUBLIC
    if(this!=&b) {free(); copy(b);}
    return *this;
  }

  _aflowlib_entry::_aflowlib_entry(const _aflowlib_entry& b) { // copy PUBLIC
    //  free();*this=b;
    copy(b);
  }

  void _aflowlib_entry::free() { // free PRIVATE
    ventry.clear();
    vauid.clear();
    vaurl.clear();
    vaflowlib_entries.clear();
    vkeywords.clear();
    vauthor.clear();
    vcorresponding.clear();
    vloop.clear();
    vcomposition.clear(); // clear all vectors
    vdft_type.clear(); // clear all vectors
    vfiles.clear(); // clear all vectors
    vfiles_LIB.clear(); // clear all vectors
    vfiles_RAW.clear(); // clear all vectors
    vfiles_WEB.clear(); // clear all vectors
    vforces.clear(); // clear all vectors
    vgeometry.clear(); // clear all vectors
    vgeometry_orig.clear(); // clear all vectors //DX20190124 - add original crystal info
    vstress_tensor.clear(); // clear all vectors
    vnbondxx.clear(); // clear all vectors
    vpositions_cartesian.clear(); // clear all vectors
    vpositions_fractional.clear(); // clear all vectors
    vspecies.clear(); // clear all vectors
    vspecies_pp.clear(); // clear all vectors
    vspecies_pp_version.clear(); // clear all vectors
    vspecies_pp_ZVAL.clear(); // clear all vectors
    vspecies_pp_AUID.clear(); // clear all vectors
    vspinD.clear(); // clear all vectors
    vspinD_magmom_orig.clear(); // clear all vectors
    vsponsor.clear();
    vstoichiometry.clear(); // clear all vectors
    vreciprocal_geometry.clear(); // clear all vectors //DX20180824 - added reciprocal lattice parameters
    // BADER
    vbader_net_charges.clear();
    vbader_atomic_volumes.clear();
    // legacy
    vserver.clear();  // clear all vectors
    for(uint i=0;i<vserverdir.size();i++)
      vserverdir.at(i).clear();
    vserverdir.clear();  // clear all vectors
    vstoich.clear(); // clear all vectors
  } 

  void _aflowlib_entry::clear() {  // clear PRIVATE
    _aflowlib_entry _temp;
    copy(_temp);
  }

  _aflowlib_entry::_aflowlib_entry(const string& file) { // constructur from file
    stringstream oss;
    if(!aurostd::FileExist(file)) {cerr << "ERROR - _aflowlib_entry::aflowlib_entr: " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " not found =" << file << endl;exit(0);} //SC20190813
    string entry;
    aurostd::efile2string(file,entry);
    Load(entry,oss);
  }

  // file2aflowlib
  uint _aflowlib_entry::file2aflowlib(const string& file,ostream& oss) {
    if(!aurostd::FileExist(file)) {cerr << "ERROR - _aflowlib_entry::file2aflowlib: " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " not found =" << file << endl;return 0;} //exit(0); //CO20170609, this is a dud
    string entry;
    aurostd::efile2string(file,entry);
    return Load(entry,oss);
  }

  // Load
  uint _aflowlib_entry::Load(const stringstream& stream,ostream& oss) {
    return Load(stream.str(),oss);
  }

  // LoadWeb
  uint _aflowlib_entry::url2aflowlib(const string& _url,ostream& oss,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "_aflowlib_entry::url2aflowlib():";
    string url=_url;
    if(url.empty()) {cerr << "ERROR - _aflowlib_entry::url2aflowlib: url.empty()" << endl;return 0;} //exit(0); //CO20170609, this is a dud
    string entry;
    if(aurostd::substring2bool(url,"index") || aurostd::substring2bool(url,"format")) {
      aurostd::StringSubst(url,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"");
      if(!aurostd::url2string(url,entry,verbose)){return 0;}   //CO, this is a dud
    } else {
      aurostd::StringSubst(url,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"");
      if(!aurostd::url2string(url+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,entry,verbose)){return 0;}  //CO, this is a dud
    }
    if(LDEBUG) {cerr << soliloquy << " entry=" << entry << endl;} //CO20180528
    return Load(entry,oss);
  }

  // Load overload
  uint _aflowlib_entry::Load(const string& _entry,ostream& oss) {
    string function = "aflowlib::_aflowlib_entry::Load()";  //ME20191119
    clear(); // start from clean
    entry=_entry; // start from loading it up !
    if(entry.empty()) {cerr << "ERROR - _aflowlib_entry::Load: entry.empty()" << endl;return 0;} //exit(0);  //CO20170609, this is a dud 
    vector<string> tokens,stokens;
    string keyword,content,line;
    aurostd::string2tokens(entry,ventry,"|");
    for(uint i=0;i<ventry.size();i++) {
      line=aurostd::RemoveWhiteSpaces(ventry.at(i));
      aurostd::string2tokens(line,tokens,"=");
      if(tokens.size()>0) {
        keyword=tokens.at(0);
        if(tokens.size()>1) {content=tokens.at(1);} else {continue;} //{content="";}  //CO20180319, content="" screws up string2double(), better to leave as AUROSTD_NAN
        if(content.empty()){continue;}  //CO20180319
        if(content=="null"){continue;}  //CO20180319 - aflux integration!
        aurostd::string2tokens(content,stokens,",");
        if(keyword=="auid") {
          auid=content; // AUID
          vauid.clear();aflowlib::auid2vauid(auid,vauid);
          // create VAUID
        }
        //CO20180409 - added the else if's for speed, no need to go through more checks than necessary
        else if(keyword=="aurl") {aurl=content;aurostd::string2tokens(content,stokens,":");for(uint j=0;j<stokens.size();j++) vaurl.push_back(stokens.at(j));}
        else if(keyword=="title") {title=content;}  //ME20190129
        else if(keyword=="keywords") {keywords=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vkeywords.push_back(stokens.at(j));}
        else if(keyword=="aflowlib_date") {aflowlib_date=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vaflowlib_date.push_back(stokens.at(j));} //CO20200624 - adding LOCK date
        else if(keyword=="aflowlib_version") {aflowlib_version=content;}
        else if(keyword=="aflowlib_entries") {aflowlib_entries=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vaflowlib_entries.push_back(stokens.at(j));}
        else if(keyword=="aflowlib_entries_number") {aflowlib_entries_number=aurostd::string2utype<int>(content);}
        else if(keyword=="aflow_version") {aflow_version=content;}
        else if(keyword=="catalog") {catalog=content;}
        else if(keyword=="data_api") {data_api=content;}
        else if(keyword=="data_source") {data_source=content;}
        else if(keyword=="data_language") {data_language=content;}
        else if(keyword=="error_status") {error_status=content;}
        else if(keyword=="author") {author=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vauthor.push_back(stokens.at(j));}
        else if(keyword=="calculation_cores") {calculation_cores=aurostd::string2utype<int>(content);}
        else if(keyword=="calculation_memory") {calculation_memory=aurostd::string2utype<double>(content);}
        else if(keyword=="calculation_time") {calculation_time=aurostd::string2utype<double>(content);}
        else if(keyword=="corresponding") {corresponding=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vcorresponding.push_back(stokens.at(j));}
        else if(keyword=="loop") {vloop.push_back(content);}  // CHECK THIS OUT IN THE FITURE
        else if(keyword=="node_CPU_Cores") {node_CPU_Cores=aurostd::string2utype<int>(content);}
        else if(keyword=="node_CPU_MHz") {node_CPU_MHz=aurostd::string2utype<double>(content);}
        else if(keyword=="node_CPU_Model") {node_CPU_Model=content;}
        else if(keyword=="node_RAM_GB") {node_RAM_GB=aurostd::string2utype<double>(content);}
        else if(keyword=="Bravais_lattice_orig") {Bravais_lattice_orig=content;}
        else if(keyword=="Bravais_lattice_relax") {Bravais_lattice_relax=content;}
        else if(keyword=="code") {code=content;}
        else if(keyword=="composition") {composition=content;for(uint j=0;j<stokens.size();j++) vcomposition.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        else if(keyword=="compound") {compound=content;}
        else if(keyword=="density") {density=aurostd::string2utype<double>(content);}
        else if(keyword=="density_orig") {density_orig=aurostd::string2utype<double>(content);} //DX20190124 - add original crystal info
        else if(keyword=="dft_type") {dft_type=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vdft_type.push_back(stokens.at(j));}
        else if(keyword=="eentropy_cell") {eentropy_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="eentropy_atom") {eentropy_atom=aurostd::string2utype<double>(content);}
        else if(keyword=="Egap") {Egap=aurostd::string2utype<double>(content);}
        else if(keyword=="Egap_fit") {Egap_fit=aurostd::string2utype<double>(content);}
        else if(keyword=="energy_cell") {energy_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="energy_atom") {energy_atom=aurostd::string2utype<double>(content);energy_atom_relax1=aurostd::string2utype<double>(content);}
        else if(keyword=="energy_cutoff") {energy_cutoff=aurostd::string2utype<double>(content);}
        else if(keyword=="delta_electronic_energy_convergence") {delta_electronic_energy_convergence=aurostd::string2utype<double>(content);}
        else if(keyword=="delta_electronic_energy_threshold") {delta_electronic_energy_threshold=aurostd::string2utype<double>(content);}
        else if(keyword=="nkpoints") {nkpoints=aurostd::string2utype<uint>(content);}
        else if(keyword=="nkpoints_irreducible") {nkpoints_irreducible=aurostd::string2utype<uint>(content);}
        else if(keyword=="kppra") {kppra=aurostd::string2utype<uint>(content);}
        else if(keyword=="kpoints") {kpoints=content;}
        else if(keyword=="kpoints_relax") {vector<int> tokens;aurostd::string2tokens(content,tokens,",");kpoints_nnn_relax=aurostd::vector2xvector(tokens);}  //ME20190129
        else if(keyword=="kpoints_static") {vector<int> tokens;aurostd::string2tokens(content,tokens,",");kpoints_nnn_static=aurostd::vector2xvector(tokens);}  //ME20190129
        else if(keyword=="kpoints_bands_path"){aurostd::string2tokens(content,kpoints_pairs,",");}  //ME20190129
        else if(keyword=="kpoints_bands_nkpts"){kpoints_bands_path_grid=aurostd::string2utype<int>(content);}  //ME20190129
        else if(keyword=="enthalpy_cell") {enthalpy_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="enthalpy_atom") {enthalpy_atom=aurostd::string2utype<double>(content);}
        else if(keyword=="enthalpy_formation_cell") {enthalpy_formation_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="enthalpy_formation_atom") {enthalpy_formation_atom=aurostd::string2utype<double>(content);}
        else if(keyword=="entropic_temperature") {entropic_temperature=aurostd::string2utype<double>(content);}
        else if(keyword=="files") {files=content;for(uint j=0;j<stokens.size();j++) vfiles.push_back(stokens.at(j));}
        else if(keyword=="files_LIB") {files_LIB=content;for(uint j=0;j<stokens.size();j++) vfiles_LIB.push_back(stokens.at(j));}
        else if(keyword=="files_RAW") {files_RAW=content;for(uint j=0;j<stokens.size();j++) vfiles_RAW.push_back(stokens.at(j));}
        else if(keyword=="files_WEB") {files_WEB=content;for(uint j=0;j<stokens.size();j++) vfiles_WEB.push_back(stokens.at(j));}
        else if(keyword=="forces") {forces=content;for(uint j=0;j<stokens.size();j++) vforces.push_back(aurostd::string2utype<double>(stokens.at(j)));}  // FIX
        else if(keyword=="geometry") {
          geometry=content; 
          vgeometry.push_back(0.0);vgeometry.push_back(0.0);vgeometry.push_back(0.0);
          vgeometry.push_back(0.0);vgeometry.push_back(0.0);vgeometry.push_back(0.0);
          if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vgeometry.at(j)=aurostd::string2utype<double>(stokens.at(j));
        }
        //DX20190124 - add original crystal info - START
        else if(keyword=="geometry_orig") {
          geometry_orig=content; 
          vgeometry_orig.push_back(0.0);vgeometry_orig.push_back(0.0);vgeometry_orig.push_back(0.0);
          vgeometry_orig.push_back(0.0);vgeometry_orig.push_back(0.0);vgeometry_orig.push_back(0.0);
          if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vgeometry_orig.at(j)=aurostd::string2utype<double>(stokens.at(j));
        }
        //DX20190124 - add original crystal info - END
        else if(keyword=="lattice_system_orig") {lattice_system_orig=content;}
        else if(keyword=="lattice_variation_orig") {lattice_variation_orig=content;}
        else if(keyword=="lattice_system_relax") {lattice_system_relax=content;}
        else if(keyword=="lattice_variation_relax") {lattice_variation_relax=content;}
        else if(keyword=="ldau_TLUJ") {ldau_TLUJ=content;}
        else if(keyword=="ldau_type") {vLDAU[0].push_back(aurostd::string2utype<double>(content));}  //ME20190129
        else if(keyword=="ldau_l") {for(uint j=0; j<stokens.size();j++) vLDAU[1].push_back(aurostd::string2utype<double>(stokens[j]));}  //ME20190129
        else if(keyword=="ldau_u") {for(uint j=0; j<stokens.size();j++) vLDAU[2].push_back(aurostd::string2utype<double>(stokens[j]));}  //ME20190129
        else if(keyword=="ldau_j") {for(uint j=0; j<stokens.size();j++) vLDAU[3].push_back(aurostd::string2utype<double>(stokens[j]));}  //ME20190129
        else if(keyword=="natoms") {natoms=aurostd::string2utype<int>(content);}
        else if(keyword=="natoms_orig") {natoms_orig=aurostd::string2utype<int>(content);} //DX20190124 - add original crystal info
        else if(keyword=="nbondxx") {nbondxx=content;for(uint j=0;j<stokens.size();j++) vnbondxx.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        else if(keyword=="nspecies") {nspecies=aurostd::string2utype<int>(content);}
        else if(keyword=="Pearson_symbol_orig") {Pearson_symbol_orig=content;}
        else if(keyword=="Pearson_symbol_relax") {Pearson_symbol_relax=content;}
        else if(keyword=="positions_cartesian") {positions_cartesian=content;for(uint j=0;j<stokens.size();j++) vpositions_cartesian.push_back(aurostd::string2utype<double>(stokens.at(j)));}  // FIX
        else if(keyword=="positions_fractional") {positions_fractional=content;for(uint j=0;j<stokens.size();j++) vpositions_fractional.push_back(aurostd::string2utype<double>(stokens.at(j)));}  // FIX
        else if(keyword=="pressure") {pressure=aurostd::string2utype<double>(content);}
        else if(keyword=="stress_tensor") {
          stress_tensor=content; 
          vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);
          vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);
          vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);
          if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vstress_tensor.at(j)=aurostd::string2utype<double>(stokens.at(j));
        }
        else if(keyword=="pressure_residual") {pressure_residual=aurostd::string2utype<double>(content);}
        else if(keyword=="Pulay_stress") {Pulay_stress=aurostd::string2utype<double>(content);}
        else if(keyword=="prototype") {prototype=content;}  // apennsy
        else if(keyword=="PV_cell") {PV_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="PV_atom") {PV_atom=aurostd::string2utype<double>(content);}
        else if(keyword=="scintillation_attenuation_length") {scintillation_attenuation_length=aurostd::string2utype<double>(content);}
        else if(keyword=="sg") {sg=content;for(uint j=0;j<stokens.size();j++) vsg.push_back(stokens.at(j));} //CO20180101
        else if(keyword=="sg2") {sg2=content;for(uint j=0;j<stokens.size();j++) vsg2.push_back(stokens.at(j));} //CO20180101
        else if(keyword=="spacegroup_orig") {spacegroup_orig=content;}
        else if(keyword=="spacegroup_relax") {spacegroup_relax=content;}
        else if(keyword=="species") {species=content;for(uint j=0;j<stokens.size();j++) vspecies.push_back(stokens.at(j));}
        else if(keyword=="species_pp") {species_pp=content;for(uint j=0;j<stokens.size();j++) vspecies_pp.push_back(stokens.at(j));}
        else if(keyword=="species_pp_version") {species_pp_version=content;for(uint j=0;j<stokens.size();j++) vspecies_pp_version.push_back(stokens.at(j));}
        else if(keyword=="species_pp_ZVAL") {species_pp_ZVAL=content;for(uint j=0;j<stokens.size();j++) vspecies_pp_ZVAL.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        else if(keyword=="species_pp_AUID") {species_pp_AUID=content;for(uint j=0;j<stokens.size();j++) vspecies_pp_AUID.push_back(stokens.at(j));}
        else if(keyword=="metagga" || keyword=="METAGGA") {METAGGA=content;}
        else if(keyword=="spin_cell") {spin_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="spin_atom") {spin_atom=aurostd::string2utype<double>(content);}
        else if(keyword=="spinD") {spinD=content;for(uint j=0;j<stokens.size();j++) vspinD.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        else if(keyword=="spinD_magmom_orig") {spinD_magmom_orig=content;for(uint j=0;j<stokens.size();j++) vspinD_magmom_orig.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        else if(keyword=="spinF") {spinF=aurostd::string2utype<double>(content);}
        else if(keyword=="sponsor") {sponsor=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vsponsor.push_back(stokens.at(j));}
        else if(keyword=="stoichiometry") {stoichiometry=content;for(uint j=0;j<stokens.size();j++) vstoichiometry.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        else if(keyword=="Egap_type") {Egap_type=content;}
        else if(keyword=="valence_cell_std") {valence_cell_std=aurostd::string2utype<double>(content);}
        else if(keyword=="valence_cell_iupac") {valence_cell_iupac=aurostd::string2utype<double>(content);}
        else if(keyword=="volume_cell") {volume_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="volume_atom") {volume_atom=aurostd::string2utype<double>(content);}
        else if(keyword=="volume_cell_orig") {volume_cell_orig=aurostd::string2utype<double>(content);} //DX20190124 - add original crystal info
        else if(keyword=="volume_atom_orig") {volume_atom_orig=aurostd::string2utype<double>(content);} //DX20190124 - add original crystal info
        // legacy
        else if(keyword=="server") {vserver.push_back(content);}
        else if(keyword=="stoich") {aurostd::string2tokens(ventry.at(i),tokens,"=");stoich=tokens.at(1);aurostd::string2tokens(stoich,stokens);for(uint j=0;j<stokens.size();j++) vstoich.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        //DX20190124 - added original symmetry info - START
        // SYMMETRY
        else if(keyword=="crystal_family_orig") {crystal_family_orig=content;}
        else if(keyword=="crystal_system_orig") {crystal_system_orig=content;}
        else if(keyword=="crystal_class_orig") {crystal_class_orig=content;}
        else if(keyword=="point_group_Hermann_Mauguin_orig") {point_group_Hermann_Mauguin_orig=content;}
        else if(keyword=="point_group_Schoenflies_orig") {point_group_Schoenflies_orig=content;}
        else if(keyword=="point_group_orbifold_orig") {point_group_orbifold_orig=content;}
        else if(keyword=="point_group_type_orig") {point_group_type_orig=content;}
        else if(keyword=="point_group_order_orig") {point_group_order=aurostd::string2utype<uint>(content);}
        else if(keyword=="point_group_structure_orig") {point_group_structure_orig=content;}
        else if(keyword=="Bravais_lattice_lattice_type_orig") {Bravais_lattice_lattice_type_orig=content;}
        else if(keyword=="Bravais_lattice_lattice_variation_type_orig") {Bravais_lattice_lattice_variation_type_orig=content;}
        else if(keyword=="Bravais_lattice_lattice_system_orig") {Bravais_lattice_lattice_system_orig=content;}
        else if(keyword=="Bravais_superlattice_lattice_type_orig") {Bravais_superlattice_lattice_type_orig=content;}
        else if(keyword=="Bravais_superlattice_lattice_variation_type_orig") {Bravais_superlattice_lattice_variation_type_orig=content;}
        else if(keyword=="Bravais_superlattice_lattice_system_orig") {Bravais_superlattice_lattice_system_orig=content;}
        else if(keyword=="Pearson_symbol_superlattice_orig") {Pearson_symbol_superlattice_orig=content;}
        else if(keyword=="reciprocal_geometry_orig") {
          reciprocal_geometry_orig=content; 
          vreciprocal_geometry_orig.push_back(0.0);vreciprocal_geometry_orig.push_back(0.0);vreciprocal_geometry_orig.push_back(0.0);
          vreciprocal_geometry_orig.push_back(0.0);vreciprocal_geometry_orig.push_back(0.0);vreciprocal_geometry_orig.push_back(0.0);
          if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vreciprocal_geometry_orig.at(j)=aurostd::string2utype<double>(stokens.at(j));
        }
        else if(keyword=="reciprocal_volume_cell_orig") {reciprocal_volume_cell_orig=aurostd::string2utype<double>(content);}
        else if(keyword=="reciprocal_lattice_type_orig") {reciprocal_lattice_type_orig=content;}
        else if(keyword=="reciprocal_lattice_variation_type_orig") {reciprocal_lattice_variation_type_orig=content;}
        else if(keyword=="Wyckoff_letters_orig") {Wyckoff_letters_orig=content;}
        else if(keyword=="Wyckoff_multiplicities_orig") {Wyckoff_multiplicities_orig=content;}
        else if(keyword=="Wyckoff_site_symmetries_orig") {Wyckoff_site_symmetries_orig=content;}
        //DX20190124 - added original symmetry info - END
        //DX20180823 - added more symmetry info - START
        // SYMMETRY
        else if(keyword=="crystal_family") {crystal_family=content;}
        else if(keyword=="crystal_system") {crystal_system=content;}
        else if(keyword=="crystal_class") {crystal_class=content;}
        else if(keyword=="point_group_Hermann_Mauguin") {point_group_Hermann_Mauguin=content;}
        else if(keyword=="point_group_Schoenflies") {point_group_Schoenflies=content;}
        else if(keyword=="point_group_orbifold") {point_group_orbifold=content;}
        else if(keyword=="point_group_type") {point_group_type=content;}
        else if(keyword=="point_group_order") {point_group_order=aurostd::string2utype<uint>(content);}
        else if(keyword=="point_group_structure") {point_group_structure=content;}
        else if(keyword=="Bravais_lattice_lattice_type") {Bravais_lattice_lattice_type=content;}
        else if(keyword=="Bravais_lattice_lattice_variation_type") {Bravais_lattice_lattice_variation_type=content;}
        else if(keyword=="Bravais_lattice_lattice_system") {Bravais_lattice_lattice_system=content;}
        else if(keyword=="Bravais_superlattice_lattice_type") {Bravais_superlattice_lattice_type=content;}
        else if(keyword=="Bravais_superlattice_lattice_variation_type") {Bravais_superlattice_lattice_variation_type=content;}
        else if(keyword=="Bravais_superlattice_lattice_system") {Bravais_superlattice_lattice_system=content;}
        else if(keyword=="Pearson_symbol_superlattice") {Pearson_symbol_superlattice=content;}
        else if(keyword=="reciprocal_geometry") {
          reciprocal_geometry=content; 
          vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);
          vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);
          if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vreciprocal_geometry.at(j)=aurostd::string2utype<double>(stokens.at(j));
        }
        else if(keyword=="reciprocal_volume_cell") {reciprocal_volume_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="reciprocal_lattice_type") {reciprocal_lattice_type=content;}
        else if(keyword=="reciprocal_lattice_variation_type") {reciprocal_lattice_variation_type=content;}
        else if(keyword=="Wyckoff_letters") {Wyckoff_letters=content;}
        else if(keyword=="Wyckoff_multiplicities") {Wyckoff_multiplicities=content;}
        else if(keyword=="Wyckoff_site_symmetries") {Wyckoff_site_symmetries=content;}
        //DX20180823 - added more symmetry info - END
        //DX20190209 - added anrl info - START
        else if(keyword=="anrl_label_orig") {anrl_label_orig=content;}
        else if(keyword=="anrl_parameter_list_orig") {anrl_parameter_list_orig=content;}
        else if(keyword=="anrl_parameter_values_orig") {anrl_parameter_values_orig=content;}
        else if(keyword=="anrl_label_relax") {anrl_label_relax=content;}
        else if(keyword=="anrl_parameter_list_relax") {anrl_parameter_list_relax=content;}
        else if(keyword=="anrl_parameter_values_relax") {anrl_parameter_values_relax=content;}
        //DX20190209 - added anrl info - END
        // AGL/AEL
        else if(keyword=="agl_thermal_conductivity_300K") {agl_thermal_conductivity_300K=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_debye") {agl_debye=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_acoustic_debye") {agl_acoustic_debye=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_gruneisen") {agl_gruneisen=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_heat_capacity_Cv_300K") {agl_heat_capacity_Cv_300K=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_heat_capacity_Cp_300K") {agl_heat_capacity_Cp_300K=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_thermal_expansion_300K") {agl_thermal_expansion_300K=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_bulk_modulus_static_300K") {agl_bulk_modulus_static_300K=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_bulk_modulus_isothermal_300K") {agl_bulk_modulus_isothermal_300K=aurostd::string2utype<double>(content);}
        else if(keyword=="agl_poisson_ratio_source") {agl_poisson_ratio_source=content;} //CT20181212
        else if(keyword=="agl_vibrational_free_energy_300K_cell") {agl_vibrational_free_energy_300K_cell=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="agl_vibrational_free_energy_300K_atom") {agl_vibrational_free_energy_300K_atom=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="agl_vibrational_entropy_300K_cell") {agl_vibrational_entropy_300K_cell=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="agl_vibrational_entropy_300K_atom") {agl_vibrational_entropy_300K_atom=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_poisson_ratio") {ael_poisson_ratio=aurostd::string2utype<double>(content);}
        else if(keyword=="ael_bulk_modulus_voigt") {ael_bulk_modulus_voigt=aurostd::string2utype<double>(content);}
        else if(keyword=="ael_bulk_modulus_reuss") {ael_bulk_modulus_reuss=aurostd::string2utype<double>(content);}
        else if(keyword=="ael_shear_modulus_voigt") {ael_shear_modulus_voigt=aurostd::string2utype<double>(content);}
        else if(keyword=="ael_shear_modulus_reuss") {ael_shear_modulus_reuss=aurostd::string2utype<double>(content);}
        else if(keyword=="ael_bulk_modulus_vrh") {ael_bulk_modulus_vrh=aurostd::string2utype<double>(content);}
        else if(keyword=="ael_shear_modulus_vrh") {ael_shear_modulus_vrh=aurostd::string2utype<double>(content);}
        else if(keyword=="ael_elastic_anisotropy") {ael_elastic_anisotropy=aurostd::string2utype<double>(content);} //CO20181129
        else if(keyword=="ael_youngs_modulus_vrh") {ael_youngs_modulus_vrh=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_speed_sound_transverse") {ael_speed_sound_transverse=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_speed_sound_longitudinal") {ael_speed_sound_longitudinal=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_speed_sound_average") {ael_speed_sound_average=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_pughs_modulus_ratio") {ael_pughs_modulus_ratio=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_debye_temperature") {ael_debye_temperature=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_applied_pressure") {ael_applied_pressure=aurostd::string2utype<double>(content);} //CT20181212
        else if(keyword=="ael_average_external_pressure") {ael_average_external_pressure=aurostd::string2utype<double>(content);} //CT20181212
        //ME20191105 BEGIN
        else if(keyword=="ael_stiffness_tensor") {
          xmatrix<double> tensor(6,6);
          vector<string> rows;
          vector<double> r;
          aurostd::string2tokens(content, rows, ";");
          if (rows.size() != 6) {
            stringstream message;
            message << "Could not read ael_stiffness_tensor: wrong number of rows"
              << " (found " << rows.size() << ", need 6).";
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
          } else {
            for (int i = 0; i < 6; i++) {
              aurostd::string2tokens(rows[i], r, ",");
              if (r.size() != 6) {
                stringstream message;
                message << "Could not read ael_stiffness_tensor: wrong number of columns"
                  << " in row " << (i + 1)
                  << " (found " << rows.size() << ", need 6).";
                throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
              } else {
                for (int j = 0; j < 6; j++) {
                  tensor[i + 1][j + 1] = r[j];
                }
              }
            }
          }
          ael_stiffness_tensor = tensor;
        } else if (keyword == "ael_compliance_tensor") {
          xmatrix<double> tensor(6,6);
          vector<string> rows;
          vector<double> r;
          aurostd::string2tokens(content, rows, ";");
          if (rows.size() != 6) {
            stringstream message;
            message << "Could not read ael_compliance_tensor: wrong number of rows"
              << " (found " << rows.size() << ", need 6).";
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
          } else {
            for (int i = 0; i < 6; i++) {
              aurostd::string2tokens(rows[i], r, ",");
              if (r.size() != 6) {
                stringstream message;
                message << "Could not read ael_compliance_tensor: wrong number of columns"
                  << " in row " << (i + 1)
                  << " (found " << rows.size() << ", need 6).";
                throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
              } else {
                for (int j = 0; j < 6; j++) {
                  tensor[i + 1][j + 1] = r[j];
                }
              }
            }
          }
          ael_compliance_tensor = tensor;
        }
        //ME20191105 END
        // BADER
        else if(keyword=="bader_net_charges") {bader_net_charges=content;aurostd::string2tokens<double>(content,vbader_net_charges,",");}
        else if(keyword=="bader_atomic_volumes") {bader_atomic_volumes=content;aurostd::string2tokens<double>(content,vbader_atomic_volumes,",");}
      }
    }
    //ME20190129 - FIX vLDAU
    if (vLDAU[0].size()) vLDAU[0].assign(vLDAU[1].size(), vLDAU[0][0]);
    // FIX LOOP
    loop="";
    vloop.push_back("aflowlib");
    //    for(uint j=0;j<vloop.size();j++) loop+=vloop.at(j)+(j<vloop.size()-1?", ":"");
    for(uint j=0;j<vloop.size();j++) loop+=vloop.at(j)+(j<vloop.size()-1?",":""); // no space
    // FIX SERVER
    server="";
    for(uint j=0;j<vserver.size();j++) {
      server+=vserver.at(j)+(j<vserver.size()-1?", ":"");
      vserverdir.push_back(*(new vector<string>(0))); // space
    }
    // FIX ICSD
    if(aurostd::substring2bool(prototype,"_ICSD_")) {
      aurostd::string2tokens(prototype,tokens,"_");
      icsd=tokens.at(tokens.size()-1);
    }
    // FIX APENNSY
    structure_name=prototype;
    structure_description=prototype;
    distance_gnd=999999; // gotta calculate it
    distance_tie=999999; // gotta calculate it
    pureA=FALSE;pureB=FALSE;
    fcc=FALSE; bcc=FALSE;hcp=FALSE;
    stoich_a=999999;stoich_b=999999;
    bond_aa=999999;bond_ab=999999;bond_bb=999999;
    vNsgroup.clear();vsgroup.clear();vstr.clear();  // apennsy
    // DONE
    if(0) {							
      bool html=TRUE;
      oss << "Keywords" << endl;
      oss << "auid=" << auid << (html?"<br>":"") << endl;
      oss << "aurl=" << aurl << (html?"<br>":"") << endl;
      oss << "title=" << title << (html?"<br>":"") << endl;  //ME20190129
      oss << "keywords=" << keywords << (html?"<br>":"") << "  vkeywords= ";for(uint j=0;j<vkeywords.size();j++) oss << vkeywords.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "Optional controls keywords (alphabetic order)" << endl;
      oss << "aflowlib_date=" << aflowlib_date << (html?"<br>":"") << endl; 
      oss << "aflowlib_version=" << aflowlib_version << (html?"<br>":"") << endl; 
      oss << "aflowlib_entries=" << aflowlib_entries << (html?"<br>":"") << endl; 
      oss << "aflowlib_entries_number=" << aflowlib_entries_number << (html?"<br>":"") << endl; 
      oss << "aflow_version=" << aflow_version << (html?"<br>":"") << endl; 
      oss << "catalog=" << catalog << (html?"<br>":"") << endl; 
      oss << "data_api=" << data_api << (html?"<br>":"") << endl; 
      oss << "data_source=" << data_source << (html?"<br>":"") << endl; 
      oss << "data_language=" << data_language << (html?"<br>":"") << endl; 
      oss << "error_status=" << error_status << (html?"<br>":"") << endl; 
      oss << "author=" << author << (html?"<br>":"") << "  vauthor= ";for(uint j=0;j<vauthor.size();j++) oss << vauthor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "calculation_cores=" << calculation_cores << (html?"<br>":"") << endl; 
      oss << "calculation_memory=" << calculation_memory << (html?"<br>":"") << endl; 
      oss << "calculation_time=" << calculation_time << (html?"<br>":"") << endl; 
      oss << "corresponding=" << corresponding << (html?"<br>":"") << "  vcorresponding= ";for(uint j=0;j<vcorresponding.size();j++) oss << vcorresponding.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "loop=" << loop << (html?"<br>":"") << "  vloop= ";for(uint j=0;j<vloop.size();j++) oss << vloop.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "node_CPU_Cores=" << node_CPU_Cores << (html?"<br>":"") << endl; 
      oss << "node_CPU_MHz=" << node_CPU_MHz << (html?"<br>":"") << endl; 
      oss << "node_CPU_Model=" << node_CPU_Model << (html?"<br>":"") << endl; 
      oss << "node_RAM_GB=" << node_RAM_GB << (html?"<br>":"") << endl; 
      oss << "Optional materials keywords (alphabetic order)" << endl;
      oss << "Bravais_lattice_orig" << Bravais_lattice_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_relax" << Bravais_lattice_relax << (html?"<br>":"") << endl;
      oss << "code=" << code << (html?"<br>":"") << endl;
      oss << "composition=" << composition << "  vcomposition= ";for(uint j=0;j<vcomposition.size();j++) oss << vcomposition.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "compound=" << compound << (html?"<br>":"") << endl;
      oss << "density=" << density << (html?"<br>":"") << endl;
      oss << "density_orig=" << density_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      oss << "dft_type=" << dft_type << (html?"<br>":"") << "  vdft_type= ";for(uint j=0;j<vdft_type.size();j++) oss << vdft_type.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "eentropy_cell=" << eentropy_cell << (html?"<br>":"") << endl; 
      oss << "eentropy_atom=" << eentropy_atom << (html?"<br>":"") << endl; 
      oss << "Egap=" << Egap << (html?"<br>":"") << endl; 
      oss << "Egap_fit=" << Egap_fit << (html?"<br>":"") << endl; 
      oss << "Egap_type=" << Egap_type << (html?"<br>":"") << endl;
      oss << "energy_cell=" << energy_cell << (html?"<br>":"") << endl; 
      oss << "energy_atom=" << energy_atom << (html?"<br>":"") << endl; 
      oss << "energy_cutoff=" << energy_cutoff << (html?"<br>":"") << endl; 
      oss << "delta_electronic_energy_convergence=" << delta_electronic_energy_convergence << (html?"<br>":"") << endl; 
      oss << "delta_electronic_energy_threshold=" << delta_electronic_energy_threshold << (html?"<br>":"") << endl; 
      oss << "nkpoints=" << nkpoints << (html?"<br>":"") << endl; 
      oss << "nkpoints_irreducible=" << nkpoints_irreducible << (html?"<br>":"") << endl; 
      oss << "kppra=" << kppra << (html?"<br>":"") << endl; 
      oss << "kpoints_relax=" << aurostd::joinWDelimiter(kpoints_nnn_relax, ",") << (html?"<br>":"") << endl;  //ME20190129
      oss << "kpoints_static=" << aurostd::joinWDelimiter(kpoints_nnn_static, ",") << (html?"<br>":"") << endl;  //ME20190129
      oss << "kpoints_bands_path=" << aurostd::joinWDelimiter(kpoints_pairs, " | ") << endl;  //ME20190129
      oss << "kpoints_bands_nkpts=" << kpoints_bands_path_grid << (html?"<br>":"") << endl;  //ME20190129
      oss << "kpoints=" << kpoints << (html?"<br>":"") << endl;      
      oss << "enthalpy_cell=" << enthalpy_cell << (html?"<br>":"") << endl; 
      oss << "enthalpy_atom=" << enthalpy_atom << (html?"<br>":"") << endl; 
      oss << "enthalpy_formation_cell=" << enthalpy_formation_cell << (html?"<br>":"") << endl; 
      oss << "enthalpy_formation_atom=" << enthalpy_formation_atom << (html?"<br>":"") << endl; 
      oss << "entropic_temperature=" << entropic_temperature << (html?"<br>":"") << endl; 
      // oss << "files=" << files << "  vfiles= ";for(uint j=0;j<vfiles.size();j++) oss << vfiles.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_LIB=" << files_LIB << "  vfiles_LIB= ";for(uint j=0;j<vfiles_LIB.size();j++) oss << vfiles_LIB.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_RAW=" << files_RAW << "  vfiles_RAW= ";for(uint j=0;j<vfiles_RAW.size();j++) oss << vfiles_RAW.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_WEB=" << files_WEB << "  vfiles_WEB= ";for(uint j=0;j<vfiles_WEB.size();j++) oss << vfiles_WEB.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "forces=" << forces << "  vforces= ";for(uint j=0;j<vforces.size();j++) oss << vforces.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "geometry=" << geometry << "  vgeometry= ";for(uint j=0;j<vgeometry.size();j++) oss << vgeometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "geometry_orig=" << geometry_orig << "  vgeometry_orig= ";for(uint j=0;j<vgeometry_orig.size();j++) oss << vgeometry_orig.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "lattice_system_orig" << lattice_system_orig << (html?"<br>":"") << endl;
      oss << "lattice_variation_orig" << lattice_variation_orig << (html?"<br>":"") << endl;
      oss << "lattice_system_relax" << lattice_system_relax << (html?"<br>":"") << endl;
      oss << "lattice_variation_relax" << lattice_variation_relax << (html?"<br>":"") << endl;
      oss << "ldau_TLUJ=" << ldau_TLUJ << (html?"<br>":"") << endl;      
      if (vLDAU[0].size()) {oss << "ldau_type=" << ((int) vLDAU[0][0]) << (html?"<br>":"") << endl;}  //ME20190129
      if (vLDAU[1].size()) {oss << "ldau_l="; oss << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[1], 0), ",") << (html?"<br>":"") << endl;}  //ME20190129
      if (vLDAU[2].size()) {oss << "ldau_u="; oss << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[2], 9), ",") << (html?"<br>":"") << endl;}  //ME20190129
      if (vLDAU[3].size()) {oss << "ldau_j="; oss << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[3], 9), ",") << (html?"<br>":"") << endl;}  //ME20190129
      oss << "natoms=" << natoms << (html?"<br>":"") << endl;
      oss << "natoms_orig=" << natoms_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      oss << "nbondxx=" << nbondxx << "  vnbondxx= ";for(uint j=0;j<vnbondxx.size();j++) oss << vnbondxx.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "nspecies=" << nspecies << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_orig" << Pearson_symbol_orig << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_relax" << Pearson_symbol_relax << (html?"<br>":"") << endl;
      // oss << "positions_cartesian=" << positions_cartesian << "  vpositions_cartesian= ";for(uint j=0;j<vpositions_cartesian.size();j++) oss << vpositions_cartesian.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "positions_fractional=" << positions_fractional << "  vpositions_fractional= ";for(uint j=0;j<vpositions_fractional.size();j++) oss << vpositions_fractional.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "pressure=" << pressure << (html?"<br>":"") << endl; 
      oss << "stress_tensor=" << stress_tensor << "  vstress_tensor= ";for(uint j=0;j<vstress_tensor.size();j++) oss << vstress_tensor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "pressure_residual=" << pressure_residual << (html?"<br>":"") << endl; 
      oss << "Pulay_stress=" << Pulay_stress << (html?"<br>":"") << endl; 
      oss << "prototype=" << prototype << (html?"<br>":"") << endl;
      oss << "PV_cell=" << PV_cell << (html?"<br>":"") << endl; 
      oss << "PV_atom=" << PV_atom << (html?"<br>":"") << endl; 
      oss << "scintillation_attenuation_length=" << scintillation_attenuation_length << (html?"<br>":"") << endl;
      oss << "sg=" << sg << (html?"<br>":"") << endl;
      oss << "sg2=" << sg2 << (html?"<br>":"") << endl;
      oss << "spacegroup_orig=" << spacegroup_orig << (html?"<br>":"") << endl;
      oss << "spacegroup_relax=" << spacegroup_relax << (html?"<br>":"") << endl;
      oss << "species=" << species << "  vspecies= ";for(uint j=0;j<vspecies.size();j++) oss << vspecies.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp=" << species_pp << "  vspecies_pp= ";for(uint j=0;j<vspecies_pp.size();j++) oss << vspecies_pp.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_version=" << species_pp_version << "  vspecies_pp_version= ";for(uint j=0;j<vspecies_pp_version.size();j++) oss << vspecies_pp_version.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_ZVAL=" << species_pp_ZVAL << "  vspecies_pp_ZVAL= ";for(uint j=0;j<vspecies_pp_ZVAL.size();j++) oss << vspecies_pp_ZVAL.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_AUID=" << species_pp_AUID << "  vspecies_pp_AUID= ";for(uint j=0;j<vspecies_pp_AUID.size();j++) oss << vspecies_pp_AUID.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "metagga=" << METAGGA << (html?"<br>":"") << endl;
      oss << "spin_cell=" << spin_cell << (html?"<br>":"") << endl; 
      oss << "spin_atom=" << spin_atom << (html?"<br>":"") << endl; 
      oss << "spinD=" << spinD << "  vspinD= "; for(uint j=0;j<vspinD.size();j++) oss << vspinD.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "spinD_magmom_orig=" << spinD_magmom_orig << "  vspinD_magmom_orig= "; for(uint j=0;j<vspinD_magmom_orig.size();j++) oss << vspinD_magmom_orig.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "spinF=" << spinF << (html?"<br>":"") << endl;
      oss << "sponsor=" << sponsor << (html?"<br>":"") << "  vsponsor= ";for(uint j=0;j<vsponsor.size();j++) oss << vsponsor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "stoichiometry=" << stoichiometry << "  vstoichiometry= ";for(uint j=0;j<vstoichiometry.size();j++) oss << vstoichiometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "valence_cell_std=" << valence_cell_std << (html?"<br>":"") << endl; 
      oss << "valence_cell_iupac=" << valence_cell_iupac << (html?"<br>":"") << endl;      
      oss << "volume_cell=" << volume_cell << (html?"<br>":"") << endl; 
      oss << "volume_atom=" << volume_atom << (html?"<br>":"") << endl; 
      oss << "volume_cell_orig=" << volume_cell_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      oss << "volume_atom_orig=" << volume_atom_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      //DX20190124 - added original symmetry info - START
      // SYMMETRY
      oss << "crystal_family_orig" << crystal_family_orig << (html?"<br>":"") << endl;
      oss << "crystal_system_orig" << crystal_system_orig << (html?"<br>":"") << endl;
      oss << "crystal_class_orig" << crystal_class_orig << (html?"<br>":"") << endl;
      oss << "point_group_Hermann_Mauguin_orig" << point_group_Hermann_Mauguin_orig << (html?"<br>":"") << endl;
      oss << "point_group_Schoenflies_orig" << point_group_Schoenflies_orig << (html?"<br>":"") << endl;
      oss << "point_group_orbifold_orig" << point_group_orbifold_orig << (html?"<br>":"") << endl;
      oss << "point_group_type_orig" << point_group_type_orig << (html?"<br>":"") << endl;
      oss << "point_group_order_orig" << point_group_order_orig << (html?"<br>":"") << endl;
      oss << "point_group_structure_orig" << point_group_structure_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_type_orig" << Bravais_lattice_lattice_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_variation_type_orig" << Bravais_lattice_lattice_variation_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_system_orig" << Bravais_lattice_lattice_system_orig << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_type_orig" << Bravais_superlattice_lattice_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_variation_type_orig" << Bravais_superlattice_lattice_variation_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_system_orig" << Bravais_superlattice_lattice_system_orig << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_superlattice_orig" << Pearson_symbol_superlattice_orig << (html?"<br>":"") << endl;
      oss << "reciprocal_geometry_orig=" << reciprocal_geometry_orig << "  vreciprocal_geometry_orig= ";for(uint j=0;j<vreciprocal_geometry_orig.size();j++) oss << vreciprocal_geometry_orig.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "reciprocal_volume_cell_orig=" << reciprocal_volume_cell_orig << (html?"<br>":"") << endl; 
      oss << "reciprocal_lattice_type_orig" << reciprocal_lattice_type_orig << (html?"<br>":"") << endl;
      oss << "reciprocal_lattice_variation_type_orig" << reciprocal_lattice_variation_type_orig << (html?"<br>":"") << endl;
      oss << "Wyckoff_letters_orig" << Wyckoff_letters_orig << (html?"<br>":"") << endl;
      oss << "Wyckoff_multiplicities_orig" << Wyckoff_multiplicities_orig << (html?"<br>":"") << endl;
      oss << "Wyckoff_site_symmetries_orig" << Wyckoff_site_symmetries_orig << (html?"<br>":"") << endl;
      //DX20190124 - added original symmetry info - END
      //DX20180823 - added more symmetry info - START
      // SYMMETRY
      oss << "crystal_family" << crystal_family << (html?"<br>":"") << endl;
      oss << "crystal_system" << crystal_system << (html?"<br>":"") << endl;
      oss << "crystal_class" << crystal_class << (html?"<br>":"") << endl;
      oss << "point_group_Hermann_Mauguin" << point_group_Hermann_Mauguin << (html?"<br>":"") << endl;
      oss << "point_group_Schoenflies" << point_group_Schoenflies << (html?"<br>":"") << endl;
      oss << "point_group_orbifold" << point_group_orbifold << (html?"<br>":"") << endl;
      oss << "point_group_type" << point_group_type << (html?"<br>":"") << endl;
      oss << "point_group_order" << point_group_order << (html?"<br>":"") << endl;
      oss << "point_group_structure" << point_group_structure << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_type" << Bravais_lattice_lattice_type << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_variation_type" << Bravais_lattice_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_system" << Bravais_lattice_lattice_system << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_type" << Bravais_superlattice_lattice_type << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_variation_type" << Bravais_superlattice_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_system" << Bravais_superlattice_lattice_system << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_superlattice" << Pearson_symbol_superlattice << (html?"<br>":"") << endl;
      oss << "reciprocal_geometry=" << reciprocal_geometry << "  vreciprocal_geometry= ";for(uint j=0;j<vreciprocal_geometry.size();j++) oss << vreciprocal_geometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "reciprocal_volume_cell=" << reciprocal_volume_cell << (html?"<br>":"") << endl; 
      oss << "reciprocal_lattice_type" << reciprocal_lattice_type << (html?"<br>":"") << endl;
      oss << "reciprocal_lattice_variation_type" << reciprocal_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Wyckoff_letters" << Wyckoff_letters << (html?"<br>":"") << endl;
      oss << "Wyckoff_multiplicities" << Wyckoff_multiplicities << (html?"<br>":"") << endl;
      oss << "Wyckoff_site_symmetries" << Wyckoff_site_symmetries << (html?"<br>":"") << endl;
      //DX20180823 - added more symmetry info - END
      //DX20190208 - added anrl info - START
      oss << "anrl_label_orig" << anrl_label_orig << (html?"<br>":"") << endl;
      oss << "anrl_parameter_list_orig" << anrl_parameter_list_orig << (html?"<br>":"") << endl;
      oss << "anrl_parameter_values_orig" << anrl_parameter_values_orig << (html?"<br>":"") << endl;
      oss << "anrl_label_relax" << anrl_label_relax << (html?"<br>":"") << endl;
      oss << "anrl_parameter_list_relax" << anrl_parameter_list_relax << (html?"<br>":"") << endl;
      oss << "anrl_parameter_values_relax" << anrl_parameter_values_relax << (html?"<br>":"") << endl;
      //DX20190208 - added anrl info - END
      // AGL/AEL
      oss << "agl_thermal_conductivity_300K" << agl_thermal_conductivity_300K << (html?"<br>":"") << endl; 
      oss << "agl_debye" << agl_debye << (html?"<br>":"") << endl; 
      oss << "agl_acoustic_debye" << agl_acoustic_debye << (html?"<br>":"") << endl; 
      oss << "agl_gruneisen" << agl_gruneisen << (html?"<br>":"") << endl; 
      oss << "agl_heat_capacity_Cv_300K" << agl_heat_capacity_Cv_300K << (html?"<br>":"") << endl; 
      oss << "agl_heat_capacity_Cp_300K" << agl_heat_capacity_Cp_300K << (html?"<br>":"") << endl; 
      oss << "agl_thermal_expansion_300K" << agl_thermal_expansion_300K << (html?"<br>":"") << endl; 
      oss << "agl_bulk_modulus_static_300K" << agl_bulk_modulus_static_300K << (html?"<br>":"") << endl; 
      oss << "agl_bulk_modulus_isothermal_300K" << agl_bulk_modulus_isothermal_300K << (html?"<br>":"") << endl; 
      oss << "agl_poisson_ratio_source" << agl_poisson_ratio_source << (html?"<br>":"") << endl; //CT20181212
      oss << "agl_vibrational_free_energy_300K_cell" << agl_vibrational_free_energy_300K_cell << (html?"<br>":"") << endl; //CT20181212
      oss << "agl_vibrational_free_energy_300K_atom" << agl_vibrational_free_energy_300K_atom << (html?"<br>":"") << endl; //CT20181212 
      oss << "agl_vibrational_entropy_300K_cell" << agl_vibrational_entropy_300K_cell << (html?"<br>":"") << endl; //CT20181212 
      oss << "agl_vibrational_entropy_300K_atom" << agl_vibrational_entropy_300K_atom << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_poisson_ratio" << ael_poisson_ratio << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_voigt" << ael_bulk_modulus_voigt << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_reuss" << ael_bulk_modulus_reuss << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_voigt" << ael_shear_modulus_voigt << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_reuss" << ael_shear_modulus_reuss << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_vrh" << ael_bulk_modulus_vrh << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_vrh" << ael_shear_modulus_vrh << (html?"<br>":"") << endl; 
      oss << "ael_elastic_anisotropy" << ael_elastic_anisotropy << (html?"<br>":"") << endl; //CO20181129
      oss << "ael_youngs_modulus_vrh" << ael_youngs_modulus_vrh << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_speed_sound_transverse" << ael_speed_sound_transverse << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_speed_sound_longitudinal" << ael_speed_sound_longitudinal << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_speed_sound_average" << ael_speed_sound_average << (html?"<br>":"") << endl; //CT20181212
      oss << "ael_pughs_modulus_ratio" << ael_pughs_modulus_ratio << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_debye_temperature" << ael_debye_temperature << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_applied_pressure" << ael_applied_pressure << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_average_external_pressure" << ael_average_external_pressure << (html?"<br>":"") << endl; //CT20181212 
      //ME20191105 BEGIN
      oss << "ael_stiffness_tensor = "; for (int i = 1; i <= 6; i++) {for (int j = 1; j <= 6; j++) oss << ael_stiffness_tensor[i][j]; oss << (html?"<br>":"") << endl;} //ME20191105
      oss << "ael_compliance_tensor = "; for (int i = 1; i <= 6; i++) {for (int j = 1; j <= 6; j++) oss << ael_compliance_tensor[i][j]; oss << (html?"<br>":"") << endl;} //ME20191105
      //ME20191105 END
      // BADER
      oss << "bader_net_charges" << bader_net_charges << "  vbader_net_charges= ";for(uint j=0;j<vbader_net_charges.size();j++) oss << vbader_net_charges.at(j) << " "; oss << (html?"<br>":"") << endl; 
      oss << "bader_atomic_volumes" << bader_atomic_volumes << "  vbader_atomic_volumes= ";for(uint j=0;j<vbader_atomic_volumes.size();j++) oss << vbader_atomic_volumes.at(j) << " "; oss << (html?"<br>":"") << endl; 
      // legacy
      oss << "server=" << server << (html?"<br>":"") << "  vserver= ";for(uint j=0;j<vserver.size();j++) oss << vserver.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "icsd=" << icsd << (html?"<br>":"") << endl;
      oss << "stoich=" << stoich << "  vstoich= ";for(uint j=0;j<vstoich.size();j++) oss << vstoich.at(j) << " "; oss << (html?"<br>":"") << endl;
    }
    return ventry.size();
  }

  // aflowlib2string 
  string _aflowlib_entry::aflowlib2string(string mode) {
    string soliloquy=XPID+"aflowlib::_aflowlib_entry::aflowlib2string():";
    stringstream sss("");
    //  string eendl="\n";

    // this is the normal aflowlib.out mode
    if(mode=="" || mode=="out" || mode=="OUT") {
      string eendl="";

      if(auid.size()) sss << "" << "aurl=" << aurl << eendl;
      if(auid.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "auid=" << auid << eendl;
      if(data_api.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_api=" << data_api << eendl;
      if(!title.empty()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "title=" << title << eendl;  //ME20190125
      if(data_source.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_source=" << data_source << eendl;
      if(data_language.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_language=" << data_language << eendl;
      if(error_status.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "error_status=" << error_status << eendl;
      // LOOP
      if(vloop.size()) {
        aurostd::sort(vloop);
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=";
        for(uint i=0;i<vloop.size();i++) sss << vloop.at(i) << (i<vloop.size()-1?",":"");
        sss << eendl;
      }
      // MATERIALS
      if(code.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "code=" << code << eendl;
      if(compound.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "compound=" << compound << eendl;
      if(prototype.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "prototype=" << prototype << eendl;
      if(nspecies!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nspecies=" << nspecies << eendl;
      if(natoms!=AUROSTD_NAN)sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "natoms=" << natoms << eendl;
      if(natoms_orig!=AUROSTD_NAN)sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "natoms_orig=" << natoms_orig << eendl; //DX20190124 - add original crystal info
      if(composition.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "composition=" << composition << eendl;
      if(density!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "density=" << density << eendl;
      if(density_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "density_orig=" << density_orig << eendl; //DX20190124 - add original crystal info
      if(scintillation_attenuation_length!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "scintillation_attenuation_length=" << scintillation_attenuation_length << eendl;
      if(stoichiometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stoichiometry=" << stoichiometry << eendl;
      if(species.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species=" << species << eendl;
      if(species_pp.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp=" << species_pp << eendl;
      if(dft_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "dft_type=" << dft_type << eendl;
      // if(species_pp_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_type=" << species_pp_type << eendl;
      if(species_pp_version.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_version=" << species_pp_version << eendl;
      if(species_pp_ZVAL.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_ZVAL=" << species_pp_ZVAL << eendl;
      if(species_pp_AUID.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_AUID=" << species_pp_AUID << eendl;
      if(METAGGA.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "metagga=" << METAGGA << eendl;
      //ME20190124 - add more detailed LDAU information
      if(ldau_TLUJ.size()) {
        //ME20190124 BEGIN
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_type=" << vLDAU[0][0] << eendl;
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_l=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[1], 0), ",") << eendl;
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_u=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[2], 9), ",") << eendl;
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_j=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[3], 9), ",") << eendl;
        //ME20190124 END
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_TLUJ=" << ldau_TLUJ << eendl;
      }
      if(valence_cell_iupac!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "valence_cell_iupac=" << valence_cell_iupac << eendl;
      if(valence_cell_std!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "valence_cell_std=" << valence_cell_std << eendl;
      if(volume_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_cell=" << volume_cell << eendl;
      if(volume_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_atom=" << volume_atom << eendl;
      if(volume_cell_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_cell_orig=" << volume_cell_orig << eendl; //DX20190124 - add original crystal info
      if(volume_atom_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_atom_orig=" << volume_atom_orig << eendl; //DX20190124 - add original crystal info
      if(pressure!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "pressure=" << pressure << eendl;
      if(stress_tensor.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stress_tensor=" << stress_tensor << eendl;
      if(pressure_residual!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "pressure_residual=" << pressure_residual << eendl;
      if(Pulay_stress!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pulay_stress=" << Pulay_stress << eendl;
      if(geometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "geometry=" << geometry << eendl;
      if(geometry_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "geometry_orig=" << geometry_orig << eendl; //DX20190124 - add original crystal info
      if(Egap!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap=" << Egap << eendl;
      if(Egap_fit!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap_fit=" << Egap_fit << eendl;
      if(Egap_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap_type=" << Egap_type << eendl;
      if(energy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_cell=" << energy_cell << eendl;
      if(energy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_atom=" << energy_atom << eendl;
      if(energy_cutoff!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_cutoff=" << energy_cutoff << eendl;
      if(delta_electronic_energy_convergence!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "delta_electronic_energy_convergence=" << delta_electronic_energy_convergence << eendl;
      if(delta_electronic_energy_threshold!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "delta_electronic_energy_threshold=" << delta_electronic_energy_threshold << eendl;
      // [NOT_PRINTED]     if(nkpoints!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nkpoints=" << nkpoints << eendl;
      // [NOT_PRINTED]     if(nkpoints_irreducible!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nkpoints_irreducible=" << nkpoints_irreducible << eendl;
      // [NOT_PRINTED]     if(kppra!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kppra=" << kppra << eendl;
      //ME20190124 BEGIN - Add the individual pieces of "kpoints" to the out file
      if ((kpoints_nnn_relax.rows == 3) && (sum(kpoints_nnn_relax) > 0)) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_relax=" << aurostd::joinWDelimiter(kpoints_nnn_relax, ",") << eendl;
      if ((kpoints_nnn_static.rows == 3) && (sum(kpoints_nnn_static) > 0)) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_static=" << aurostd::joinWDelimiter(kpoints_nnn_static, ",") << eendl;
      if (kpoints_pairs.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_bands_path=" << aurostd::joinWDelimiter(kpoints_pairs, ",") << eendl;
      if (kpoints_bands_path_grid > 0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_bands_nkpts=" << ((int) kpoints_bands_path_grid) << eendl;
      //ME20190124 END
      if(kpoints.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints=" << kpoints << eendl;
      if(enthalpy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_cell=" << enthalpy_cell << eendl;
      if(enthalpy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_atom=" << enthalpy_atom << eendl;
      if(eentropy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "eentropy_cell=" << eentropy_cell << eendl;
      if(eentropy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "eentropy_atom=" << eentropy_atom << eendl;
      if(enthalpy_formation_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_cell=" << enthalpy_formation_cell << eendl;
      if(enthalpy_formation_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_atom=" << enthalpy_formation_atom << eendl;
      if(entropic_temperature!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "entropic_temperature=" << entropic_temperature << eendl;
      if(PV_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "PV_cell=" << PV_cell << eendl;
      if(PV_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "PV_atom=" << PV_atom << eendl;
      if(spin_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spin_cell=" << spin_cell << eendl;
      if(spin_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spin_atom=" << spin_atom << eendl;
      if(spinD.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinD=" << spinD << eendl;
      if(spinD_magmom_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinD_magmom_orig=" << spinD_magmom_orig << eendl;
      if(spinF!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinF=" << spinF << eendl;
      if(stoich.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stoich=" << stoich << eendl;
      if(calculation_time!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_time=" << calculation_time << eendl;
      if(calculation_memory!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_memory=" << calculation_memory << eendl;
      if(calculation_cores!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_cores=" << calculation_cores << eendl;
      if(nbondxx.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nbondxx=" << nbondxx << eendl;
      if(sg.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "sg=" << sg << eendl;
      if(sg2.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "sg2=" << sg2 << eendl;
      if(spacegroup_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spacegroup_orig=" << spacegroup_orig << eendl;
      if(spacegroup_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spacegroup_relax=" << spacegroup_relax << eendl;
      if(forces.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "forces=" << forces << eendl;
      if(positions_cartesian.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "positions_cartesian=" << positions_cartesian << eendl;
      if(positions_fractional.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "positions_fractional=" << positions_fractional << eendl;
      if(Bravais_lattice_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_orig=" << Bravais_lattice_orig << eendl;
      if(lattice_variation_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_variation_orig=" << lattice_variation_orig << eendl;
      if(lattice_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_system_orig=" << lattice_system_orig << eendl;
      if(Pearson_symbol_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_orig=" << Pearson_symbol_orig << eendl;
      if(Bravais_lattice_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_relax=" << Bravais_lattice_relax << eendl;
      if(lattice_variation_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_variation_relax=" << lattice_variation_relax << eendl;
      if(lattice_system_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_system_relax=" << lattice_system_relax << eendl;
      if(Pearson_symbol_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_relax=" << Pearson_symbol_relax << eendl;
      //DX20190124 - added original symmetry info - START
      // SYMMETRY
      if(crystal_family_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_family_orig=" << crystal_family_orig << eendl;
      if(crystal_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_system_orig=" << crystal_system_orig << eendl;
      if(crystal_class_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_class_orig=" << crystal_class_orig << eendl;
      if(point_group_Hermann_Mauguin_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Hermann_Mauguin_orig=" << point_group_Hermann_Mauguin_orig << eendl;
      if(point_group_Schoenflies_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Schoenflies_orig=" << point_group_Schoenflies_orig << eendl;
      if(point_group_orbifold_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_orbifold_orig=" << point_group_orbifold_orig << eendl;
      if(point_group_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_type_orig=" << point_group_type_orig << eendl;
      if(point_group_order_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_order_orig=" << point_group_order_orig << eendl;
      if(point_group_structure_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_structure_orig=" << point_group_structure_orig << eendl;
      if(Bravais_lattice_lattice_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_type_orig=" << Bravais_lattice_lattice_type_orig << eendl;
      if(Bravais_lattice_lattice_variation_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_variation_type_orig=" << Bravais_lattice_lattice_variation_type_orig << eendl;
      if(Bravais_lattice_lattice_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_system_orig=" << Bravais_lattice_lattice_system_orig << eendl;
      if(Bravais_superlattice_lattice_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_type_orig=" << Bravais_superlattice_lattice_type_orig << eendl;
      if(Bravais_superlattice_lattice_variation_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_variation_type_orig=" << Bravais_superlattice_lattice_variation_type_orig << eendl;
      if(Bravais_superlattice_lattice_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_system_orig=" << Bravais_superlattice_lattice_system_orig << eendl;
      if(Pearson_symbol_superlattice_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_superlattice_orig=" << Pearson_symbol_superlattice_orig << eendl;
      if(reciprocal_geometry_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_geometry_orig=" << reciprocal_geometry_orig << eendl;
      if(reciprocal_volume_cell_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_volume_cell_orig=" << reciprocal_volume_cell_orig << eendl;
      if(reciprocal_lattice_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_type_orig=" << reciprocal_lattice_type_orig << eendl;
      if(reciprocal_lattice_variation_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_variation_type_orig=" << reciprocal_lattice_variation_type_orig << eendl;
      if(Wyckoff_letters_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_letters_orig=" << Wyckoff_letters_orig << eendl;
      if(Wyckoff_multiplicities_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_multiplicities_orig=" << Wyckoff_multiplicities_orig << eendl;
      if(Wyckoff_site_symmetries_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_site_symmetries_orig=" << Wyckoff_site_symmetries_orig << eendl;
      //DX20190124 - added original symmetry info - END
      //DX20180823 - added more symmetry info - START
      // SYMMETRY
      if(crystal_family.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_family=" << crystal_family << eendl;
      if(crystal_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_system=" << crystal_system << eendl;
      if(crystal_class.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_class=" << crystal_class << eendl;
      if(point_group_Hermann_Mauguin.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Hermann_Mauguin=" << point_group_Hermann_Mauguin << eendl;
      if(point_group_Schoenflies.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Schoenflies=" << point_group_Schoenflies << eendl;
      if(point_group_orbifold.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_orbifold=" << point_group_orbifold << eendl;
      if(point_group_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_type=" << point_group_type << eendl;
      if(point_group_order!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_order=" << point_group_order << eendl;
      if(point_group_structure.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_structure=" << point_group_structure << eendl;
      if(Bravais_lattice_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_type=" << Bravais_lattice_lattice_type << eendl;
      if(Bravais_lattice_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_variation_type=" << Bravais_lattice_lattice_variation_type << eendl;
      if(Bravais_lattice_lattice_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_system=" << Bravais_lattice_lattice_system << eendl;
      if(Bravais_superlattice_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_type=" << Bravais_superlattice_lattice_type << eendl;
      if(Bravais_superlattice_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_variation_type=" << Bravais_superlattice_lattice_variation_type << eendl;
      if(Bravais_superlattice_lattice_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_system=" << Bravais_superlattice_lattice_system << eendl;
      if(Pearson_symbol_superlattice.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_superlattice=" << Pearson_symbol_superlattice << eendl;
      if(reciprocal_geometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_geometry=" << reciprocal_geometry << eendl;
      if(reciprocal_volume_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_volume_cell=" << reciprocal_volume_cell << eendl;
      if(reciprocal_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_type=" << reciprocal_lattice_type << eendl;
      if(reciprocal_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_variation_type=" << reciprocal_lattice_variation_type << eendl;
      if(Wyckoff_letters.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_letters=" << Wyckoff_letters << eendl;
      if(Wyckoff_multiplicities.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_multiplicities=" << Wyckoff_multiplicities << eendl;
      if(Wyckoff_site_symmetries.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_site_symmetries=" << Wyckoff_site_symmetries << eendl;
      //DX20180823 - added more symmetry info - END
      //DX20190208 - added anrl info - START
      if(anrl_label_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "anrl_label_orig=" << anrl_label_orig << eendl;
      if(anrl_parameter_list_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "anrl_parameter_list_orig=" << anrl_parameter_list_orig << eendl;
      if(anrl_parameter_values_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "anrl_parameter_values_orig=" << anrl_parameter_values_orig << eendl;
      if(anrl_label_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "anrl_label_relax=" << anrl_label_relax << eendl;
      if(anrl_parameter_list_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "anrl_parameter_list_relax=" << anrl_parameter_list_relax << eendl;
      if(anrl_parameter_values_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "anrl_parameter_values_relax=" << anrl_parameter_values_relax << eendl;
      //DX20190208 - added anrl info - END
      // AGL/AEL
      if(agl_thermal_conductivity_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_thermal_conductivity_300K=" << agl_thermal_conductivity_300K << eendl;
      if(agl_debye!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_debye=" << agl_debye << eendl;
      if(agl_acoustic_debye!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_acoustic_debye=" << agl_acoustic_debye << eendl;
      if(agl_gruneisen!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_gruneisen=" << agl_gruneisen << eendl;
      if(agl_heat_capacity_Cv_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_heat_capacity_Cv_300K=" << agl_heat_capacity_Cv_300K << eendl;
      if(agl_heat_capacity_Cp_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_heat_capacity_Cp_300K=" << agl_heat_capacity_Cp_300K << eendl;
      if(agl_thermal_expansion_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_thermal_expansion_300K=" << agl_thermal_expansion_300K << eendl;
      if(agl_bulk_modulus_static_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_bulk_modulus_static_300K=" << agl_bulk_modulus_static_300K << eendl;
      if(agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_bulk_modulus_isothermal_300K=" << agl_bulk_modulus_isothermal_300K << eendl;
      if(agl_poisson_ratio_source.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_poisson_ratio_source=" << agl_poisson_ratio_source << eendl; //CT20181212
      if(agl_vibrational_free_energy_300K_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_free_energy_300K_cell=" << agl_vibrational_free_energy_300K_cell << eendl; //CT20181212
      if(agl_vibrational_free_energy_300K_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_free_energy_300K_atom=" << agl_vibrational_free_energy_300K_atom << eendl; //CT20181212
      if(agl_vibrational_entropy_300K_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_entropy_300K_cell=" << agl_vibrational_entropy_300K_cell << eendl; //CT20181212
      if(agl_vibrational_entropy_300K_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_entropy_300K_atom=" << agl_vibrational_entropy_300K_atom << eendl; //CT20181212
      if(ael_poisson_ratio!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_poisson_ratio=" << ael_poisson_ratio << eendl;
      if(ael_bulk_modulus_voigt!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_voigt=" << ael_bulk_modulus_voigt << eendl;
      if(ael_bulk_modulus_reuss!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_reuss=" << ael_bulk_modulus_reuss << eendl;
      if(ael_shear_modulus_voigt!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_voigt=" << ael_shear_modulus_voigt << eendl;
      if(ael_shear_modulus_reuss!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_reuss=" << ael_shear_modulus_reuss << eendl;
      if(ael_bulk_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_vrh=" << ael_bulk_modulus_vrh << eendl;
      if(ael_shear_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_vrh=" << ael_shear_modulus_vrh << eendl;
      if(ael_elastic_anisotropy!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_elastic_anisotropy=" << ael_elastic_anisotropy << eendl; //CO20181129
      if(ael_youngs_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_youngs_modulus_vrh=" << ael_youngs_modulus_vrh << eendl; //CT20181212
      if(ael_speed_sound_transverse!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_speed_sound_transverse=" << ael_speed_sound_transverse << eendl; //CT20181212
      if(ael_speed_sound_longitudinal!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_speed_sound_longitudinal=" << ael_speed_sound_longitudinal << eendl; //CT20181212
      if(ael_speed_sound_average!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_speed_sound_average=" << ael_speed_sound_average << eendl; //CT20181212
      if(ael_pughs_modulus_ratio!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_pughs_modulus_ratio=" << ael_pughs_modulus_ratio << eendl; //CT20181212
      if(ael_debye_temperature!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_debye_temperature=" << ael_debye_temperature << eendl; //CT20181212
      if(ael_applied_pressure!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_applied_pressure=" << ael_applied_pressure << eendl; //CT20181212
      if(ael_average_external_pressure!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_average_external_pressure=" << ael_average_external_pressure << eendl; //CT20181212
      //ME20191105 BEGIN
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_stiffness_tensor=";
      if ((ael_stiffness_tensor.rows == 6) && (ael_stiffness_tensor.cols == 6)) {
        for (int i = 1; i <= 6; i++) {
          for (int j = 1; j <= 6; j++) sss << ael_stiffness_tensor[i][j] << ((j < 6)?",":"");
          sss << ((i < 6)?";":"") << eendl;
        }
      }
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_compliance_tensor=";
      if ((ael_compliance_tensor.rows == 6) && (ael_compliance_tensor.cols == 6)) {
        for (int i = 1; i <= 6; i++) {
          for (int j = 1; j <= 6; j++) {
            sss << ael_compliance_tensor[i][j] << ((j < 6)?",":"");
          }
          sss << ((i < 6)?";":"") << eendl;
        }
      }
      //ME20191105 END
      // BADER
      if(bader_net_charges.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "bader_net_charges=" << bader_net_charges << eendl;
      if(bader_atomic_volumes.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "bader_atomic_volumes=" << bader_atomic_volumes << eendl;
      // FILES
      if(vfiles.size()) {
        aurostd::sort(vfiles);
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "files=";
        for(uint i=0;i<vfiles.size();i++) sss << vfiles.at(i) << (i<vfiles.size()-1?",":"");
        sss << eendl;
      }
      // CPUS
      if(node_CPU_Model.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_Model=" << node_CPU_Model << eendl;
      if(node_CPU_Cores!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_Cores=" << node_CPU_Cores << eendl;
      if(node_CPU_MHz!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_MHz=" << node_CPU_MHz << eendl;
      if(node_RAM_GB!=INF) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_RAM_GB=" << node_RAM_GB << eendl;
      // VERSION/DATE
      if(aflow_version.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_version=" << aflow_version << eendl;
      if(catalog.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "catalog=" << catalog << eendl;
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflowlib_version=" << string(AFLOW_VERSION) << eendl;
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflowlib_date=" << aurostd::joinWDelimiter(vaflowlib_date,",") << eendl;  //CO20200624 - adding LOCK date
      sss << endl;

    } // out

    // this is the aflowlib.json mode
    if(mode=="json" || mode=="JSON") {  //CO OPERATE HERE ALL THE STRINGS AS BEFORE
      string eendl="";
      bool PRINT_NULL=FALSE;
      stringstream sscontent_json;
      vector<string> vcontent_json;
      vector<string> sg_tokens;
      stringstream ss_helper;
      vector<vector<string> > vvs;
      vector<string> vs;
      bool odd_xvec_count;

      //////////////////////////////////////////////////////////////////////////
      if(auid.size()) {
        sscontent_json << "\"aurl\":\"" << aurl << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"aurl\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(auid.size()) {
        sscontent_json << "\"auid\":\"" << auid << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"auid\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //ME20190125 BEGIN
      //////////////////////////////////////////////////////////////////////////
      if(!title.empty()) {
        sscontent_json << "\"title\":\"" << title << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"title\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      //ME20190125 END

      //////////////////////////////////////////////////////////////////////////
      if(data_api.size()) {
        sscontent_json << "\"data_api\":\"" << data_api << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_api\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(data_source.size()) {
        sscontent_json << "\"data_source\":\"" << data_source << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_source\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(data_language.size()) {
        sscontent_json << "\"data_language\":\"" << data_language << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_language\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(error_status.size()) {
        sscontent_json << "\"error_status\":\"" << error_status << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"error_status\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // LOOP
      //////////////////////////////////////////////////////////////////////////
      if(vloop.size()) {
        aurostd::sort(vloop);
        sscontent_json << "\"loop\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vloop,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"loop\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // MATERIALS
      //////////////////////////////////////////////////////////////////////////
      if(code.size()) {
        sscontent_json << "\"code\":\"" << code << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"code\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(compound.size()) {
        sscontent_json << "\"compound\":\"" << compound << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"compound\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(prototype.size()) {
        sscontent_json << "\"prototype\":\"" << prototype << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"prototype\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(nspecies!=AUROSTD_NAN) {
        sscontent_json << "\"nspecies\":" << nspecies << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"nspecies\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(natoms!=AUROSTD_NAN) {
        sscontent_json << "\"natoms\":" << natoms << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"natoms\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(natoms_orig!=AUROSTD_NAN) {
        sscontent_json << "\"natoms_orig\":" << natoms_orig << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"natoms_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(vcomposition.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"composition\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vcomposition,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"composition\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(density!=AUROSTD_NAN) {
        sscontent_json << "\"density\":" << density << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"density\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(density_orig!=AUROSTD_NAN) {
        sscontent_json << "\"density_orig\":" << density_orig << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"density_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(scintillation_attenuation_length!=AUROSTD_NAN) {
        sscontent_json << "\"scintillation_attenuation_length\":" << scintillation_attenuation_length << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"scintillation_attenuation_length\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vstoichiometry.size()) {
        //aflowlib_libraries specifies precision of 9
        sscontent_json << "\"stoichiometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstoichiometry,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"stoichiometry\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies.size()) {
        sscontent_json << "\"species\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp.size()) {
        sscontent_json << "\"species_pp\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(dft_type.size()) {
        //DX+CO START
        //sscontent_json << "\"dft_type\":\"" << dft_type << "\"" << eendl; //CO, this is technically a vector (RESTAPI paper)
        sscontent_json << "\"dft_type\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vdft_type,"\""),",") << "]" << eendl;
        //DX+CO END
      } else {
        if(PRINT_NULL) sscontent_json << "\"dft_type\":null" << dft_type << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // if(species_pp_type.size()) sscontent_json << "species_pp_type=" << species_pp_type << eendl;

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_version.size()) {
        sscontent_json << "\"species_pp_version\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp_version,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_version\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_ZVAL.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"species_pp_ZVAL\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspecies_pp_ZVAL,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_ZVAL\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_AUID.size()) {
        sscontent_json << "\"species_pp_AUID\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp_AUID,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_AUID\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(METAGGA.size()) {
        sscontent_json << "\"metagga\":\"" << METAGGA << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"metagga\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //ME20190124 - Modified to include more detailed LDAU information
      if(ldau_TLUJ.size()) {
        ss_helper.str("");
        vs.clear();
        //only string is available, so we have to parse really fast
        //would be nice if we have vldau_TLUJ already
        //ME20190124 - vLDAU is now available
        //[OBSOLETE ME20190124 - use vLDAU] vector<string> ldau_TLUJ_tokens;
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ,ldau_TLUJ_tokens,";");
        //[OBSOLETE ME20190124 - use vLDAU] if(ldau_TLUJ_tokens.size()==4){
        //[OBSOLETE ME20190124 - use vLDAU] conversion to double ENSURES that these are numbers
        //[OBSOLETE ME20190124 - use vLDAU] non-numbers without "" will break json
        //[OBSOLETE ME20190124 - use vLDAU] int T=aurostd::string2utype<int>(ldau_TLUJ_tokens.at(0));
        int T = (int) vLDAU[0][0]; //ME20190124
        vector<int> L(vLDAU[1].begin(), vLDAU[1].end());  //ME20190124
        vector<double> U = vLDAU[2],J = vLDAU[3];  //ME20190124
        //[OBSOLETE ME20190124 - NOT USED] vector<string> ldau_TLUJ_tokens2;
        //[OBSOLETE ME20190124 - use vLDAU] breaking up not necessary, but a nice check that we don't have hanging commas
        //[OBSOLETE ME20190124 - use vLDAU] the extra space at the end will be removed by joinWDelimiter()
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ_tokens.at(1),L,",");
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ_tokens.at(2),U,",");
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ_tokens.at(3),J,",");
        if(L.size()&&U.size()&&J.size()){
          //no precision needed
          vs.push_back(aurostd::utype2string(T));
          vs.push_back("["+aurostd::joinWDelimiter(L,",")+"]");
          vs.push_back("["+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(U,9),",")+"]");
          vs.push_back("["+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(J,9),",")+"]");
          ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
          vector<string> ldau_keys;  //ME20190124
          aurostd::string2tokens("ldau_type,ldau_l,ldau_u,ldau_j", ldau_keys, ",");  //ME20190124
          for (uint i = 0; i < ldau_keys.size(); i++) {  //ME20190124
            sscontent_json << "\"" << ldau_keys[i] << "\":" << vs[i] << eendl;  //ME20190124
            vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");  //ME20190124
          }
          vs.clear();
        }
        //} ME20190124
        vs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"ldau_TLUJ\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"ldau_TLUJ\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(valence_cell_iupac!=AUROSTD_NAN) {
        sscontent_json << "\"valence_cell_iupac\":" << valence_cell_iupac << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"valence_cell_iupac\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(valence_cell_std!=AUROSTD_NAN) {
        sscontent_json << "\"valence_cell_std\":" << valence_cell_std << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"valence_cell_std\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_cell!=AUROSTD_NAN) {
        sscontent_json << "\"volume_cell\":" << volume_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_atom!=AUROSTD_NAN) {
        sscontent_json << "\"volume_atom\":" << volume_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(volume_cell_orig!=AUROSTD_NAN) {
        sscontent_json << "\"volume_cell_orig\":" << volume_cell_orig << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_cell_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_atom_orig!=AUROSTD_NAN) {
        sscontent_json << "\"volume_atom_orig\":" << volume_atom_orig << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_atom_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(pressure!=AUROSTD_NAN) {
        sscontent_json << "\"pressure\":" << pressure << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"pressure\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vstress_tensor.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"stress_tensor\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstress_tensor,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"stress_tensor\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(pressure_residual!=AUROSTD_NAN) {
        sscontent_json << "\"pressure_residual\":" << pressure_residual << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"pressure_residual\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pulay_stress!=AUROSTD_NAN) {
        sscontent_json << "\"Pulay_stress\":" << Pulay_stress << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pulay_stress\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vgeometry.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"geometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vgeometry,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"geometry\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(vgeometry_orig.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"geometry_orig\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vgeometry_orig,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"geometry_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(Egap!=AUROSTD_NAN) {
        sscontent_json << "\"Egap\":" << Egap << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Egap_fit!=AUROSTD_NAN) {
        sscontent_json << "\"Egap_fit\":" << Egap_fit << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap_fit\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Egap_type.size()) {
        sscontent_json << "\"Egap_type\":\"" << Egap_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"energy_cell\":" << energy_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"energy_atom\":" << energy_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_cutoff!=AUROSTD_NAN) {
        sscontent_json << "\"energy_cutoff\":" << energy_cutoff << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_cutoff\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(delta_electronic_energy_convergence!=AUROSTD_NAN) {
        sscontent_json << "\"delta_electronic_energy_convergence\":" << delta_electronic_energy_convergence << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"delta_electronic_energy_convergence\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(delta_electronic_energy_threshold!=AUROSTD_NAN) {
        sscontent_json << "\"delta_electronic_energy_threshold\":" << delta_electronic_energy_threshold << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"delta_electronic_energy_threshold\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(nkpoints!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"nkpoints\":" << nkpoints << eendl;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"nkpoints\":null" << eendl;
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(nkpoints_irreducible!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"nkpoints_irreducible\":" << nkpoints_irreducible << eendl;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"nkpoints_irreducible\":null" << eendl;
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(kppra!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"kppra\":" << kppra << eendl;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"kppra\":null" << eendl;
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(kpoints.size()) {
        //this one is a bit complicated, so we will test if the string was created, and then recreate the json array on the spot
        ss_helper.str("");
        vs.clear();
        //ME20190124 - Add the individual pieces of "kpoints" to the json file
        if ((kpoints_nnn_relax.rows==3) && (sum(kpoints_nnn_relax) > 0)) {  //ME20190128
          vs.push_back("["+aurostd::joinWDelimiter(kpoints_nnn_relax,",")+"]");
          sscontent_json << "\"kpoints_relax\":" << vs.back() << eendl;  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_relax\":null" << eendl;  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");  //ME20190124
        if ((kpoints_nnn_static.rows==3) && (sum(kpoints_nnn_static) > 0)) {  //ME20190128
          vs.push_back("["+aurostd::joinWDelimiter(kpoints_nnn_static,",")+"]");
          if (sum(kpoints_nnn_static) > 0) sscontent_json << "\"kpoints_static\":" << vs.back() << eendl;  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_static\":null" << eendl;  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");  //ME20190124
        if(kpoints_pairs.size()){
          //first for escape characters in \Gamma or \Sigma
          vector<string> kpoints_pairs_new;
          char issue_c='\\';
          stringstream issue_ss; issue_ss << issue_c;
          string fix_s="\\\\";
          for(uint i=0;i<kpoints_pairs.size();i++){
            kpoints_pairs_new.push_back(aurostd::StringSubst(kpoints_pairs.at(i),issue_ss.str(),fix_s));
          }
          vs.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(kpoints_pairs_new,"\""),",")+"]");
          sscontent_json << "\"kpoints_bands_path\":" << vs.back() << eendl;  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_bands_path\":null" << eendl;  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");  //ME20190124
        ss_helper << aurostd::joinWDelimiter(vs, ",");  //ME20190128
        if(kpoints_bands_path_grid!=0){
          //ME20190128 - This causes kpoints to only be written when the band structure
          // was calculated. This is inconsistent with the aflowlib.out file
          // [OBSOLETE ME20190128] ss_helper << aurostd::joinWDelimiter(vs,",") << "," << aurostd::utype2string(kpoints_bands_path_grid);
          ss_helper << "," << aurostd::utype2string(kpoints_bands_path_grid);
          sscontent_json << "\"kpoints_bands_nkpts\":" << ((int) kpoints_bands_path_grid) << eendl;  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_bands_nkpts\":null" << eendl;  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");  //ME20190124
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"kpoints\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"kpoints\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_cell\":" << enthalpy_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_atom\":" << enthalpy_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(eentropy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"eentropy_cell\":" << eentropy_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"eentropy_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(eentropy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"eentropy_atom\":" << eentropy_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"eentropy_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_cell!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_formation_cell\":" << enthalpy_formation_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_atom!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_formation_atom\":" << enthalpy_formation_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(entropic_temperature!=AUROSTD_NAN) {
        sscontent_json << "\"entropic_temperature\":" << entropic_temperature << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"entropic_temperature\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(PV_cell!=AUROSTD_NAN) {
        sscontent_json << "\"PV_cell\":" << PV_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"PV_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(PV_atom!=AUROSTD_NAN) {
        sscontent_json << "\"PV_atom\":" << PV_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"PV_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spin_cell!=AUROSTD_NAN) {
        sscontent_json << "\"spin_cell\":" << spin_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spin_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spin_atom!=AUROSTD_NAN) {
        sscontent_json << "\"spin_atom\":" << spin_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spin_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspinD.size()) {
        //aflowlib_libraries specifies precision of 5
        sscontent_json << "\"spinD\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspinD,5),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinD\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspinD_magmom_orig.size()) {
        //aflowlib_libraries specifies precision of 5
        sscontent_json << "\"spinD_magmom_orig\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspinD_magmom_orig,5),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinD_magmom_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spinF!=AUROSTD_NAN) {
        sscontent_json << "\"spinF\":" << spinF << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinF\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // [OBSOLETE] 
      // [OBSOLETE] if(stoich.size()) {
      // [OBSOLETE]   //just use the string SC made
      // [OBSOLETE]   sscontent_json << "\"stoich\":\"" << stoich << "\"" << eendl;
      // [OBSOLETE]   ////aflowlib_libraries does not specify precision
      // [OBSOLETE]   //sscontent_json << "\"stoich\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstoich,9),",") << "]" << eendl;
      // [OBSOLETE] } else {
      // [OBSOLETE]   if(PRINT_NULL) sscontent_json << "\"stoich\":null" << eendl;
      // [OBSOLETE] }
      // [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_time!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_time\":" << calculation_time << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_time\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_memory!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_memory\":" << calculation_memory << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_memory\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_cores!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_cores\":" << calculation_cores << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_cores\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vnbondxx.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"nbondxx\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vnbondxx,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"nbondxx\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(sg.size()) {
        aurostd::string2tokens(sg,sg_tokens,",");
        sscontent_json << "\"sg\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(sg_tokens,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"sg\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(sg2.size()) {
        aurostd::string2tokens(sg2,sg_tokens,",");
        sscontent_json << "\"sg2\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(sg_tokens,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"sg2\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spacegroup_orig.size()) {
        sscontent_json << "\"spacegroup_orig\":" << spacegroup_orig << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spacegroup_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spacegroup_relax.size()) {
        sscontent_json << "\"spacegroup_relax\":" << spacegroup_relax << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spacegroup_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vforces.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vforces.size();i++){
          if(vforces.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vforces.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"forces\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"forces\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vpositions_cartesian.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vpositions_cartesian.size();i++){
          if(vpositions_cartesian.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vpositions_cartesian.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"positions_cartesian\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"positions_cartesian\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vpositions_fractional.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vpositions_fractional.size();i++){
          if(vpositions_fractional.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vpositions_fractional.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"positions_fractional\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"positions_fractional\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_orig.size()) {
        sscontent_json << "\"Bravais_lattice_orig\":\"" << Bravais_lattice_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_variation_orig.size()) {
        sscontent_json << "\"lattice_variation_orig\":\"" << lattice_variation_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_variation_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_system_orig.size()) {
        sscontent_json << "\"lattice_system_orig\":\"" << lattice_system_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_system_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_orig.size()) {
        sscontent_json << "\"Pearson_symbol_orig\":\"" << Pearson_symbol_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_relax.size()) {
        sscontent_json << "\"Bravais_lattice_relax\":\"" << Bravais_lattice_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_variation_relax.size()) {
        sscontent_json << "\"lattice_variation_relax\":\"" << lattice_variation_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_variation_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_system_relax.size()) {
        sscontent_json << "\"lattice_system_relax\":\"" << lattice_system_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_system_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_relax.size()) {
        sscontent_json << "\"Pearson_symbol_relax\":\"" << Pearson_symbol_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - added original symmetry info - START
      // SYMMETRY
      //////////////////////////////////////////////////////////////////////////
      if(crystal_family_orig.size()){
        sscontent_json << "\"crystal_family_orig\":\"" << crystal_family_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_family_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(crystal_system_orig.size()){
        sscontent_json << "\"crystal_system_orig\":\"" << crystal_system_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_system_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(crystal_class_orig.size()){
        sscontent_json << "\"crystal_class_orig\":\"" << crystal_class_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL){ sscontent_json << "\"crystal_class_orig\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Hermann_Mauguin_orig.size()){
        sscontent_json << "\"point_group_Hermann_Mauguin_orig\":\"" << point_group_Hermann_Mauguin_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Hermann_Mauguin_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Schoenflies_orig.size()){
        sscontent_json << "\"point_group_Schoenflies_orig\":\"" << point_group_Schoenflies_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Schoenflies_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_orbifold_orig.size()){
        sscontent_json << "\"point_group_orbifold_orig\":\"" << point_group_orbifold_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_orbifold_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_type_orig.size()){
        sscontent_json << "\"point_group_type_orig\":\"" << point_group_type_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_type_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_order_orig!=AUROSTD_NAN){
        sscontent_json << "\"point_group_order_orig\":" << point_group_order_orig << eendl; //DX20190124 - changed to number, not string 
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_order_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_structure_orig.size()){
        sscontent_json << "\"point_group_structure_orig\":\"" << point_group_structure_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_structure_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_type_orig.size()){
        sscontent_json << "\"Bravais_lattice_lattice_type_orig\":\"" << Bravais_lattice_lattice_type_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_type_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_variation_type_orig.size()){
        sscontent_json << "\"Bravais_lattice_lattice_variation_type_orig\":\"" << Bravais_lattice_lattice_variation_type_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_variation_type_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_system_orig.size()){
        sscontent_json << "\"Bravais_lattice_lattice_system_orig\":\"" << Bravais_lattice_lattice_system_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_system_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_type_orig.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_type_orig\":\"" << Bravais_superlattice_lattice_type_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_type_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_variation_type_orig.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_variation_type_orig\":\"" << Bravais_superlattice_lattice_variation_type_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_variation_type_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_system_orig.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_system_orig\":\"" << Bravais_superlattice_lattice_system_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_system_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_superlattice_orig.size()){
        sscontent_json << "\"Pearson_symbol_superlattice_orig\":\"" << Pearson_symbol_superlattice_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_superlattice_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(vreciprocal_geometry_orig.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"reciprocal_geometry_orig\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vreciprocal_geometry_orig,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_geometry_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_volume_cell_orig!=AUROSTD_NAN) {
        sscontent_json << "\"reciprocal_volume_cell_orig\":" << reciprocal_volume_cell_orig << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_volume_cell_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_type_orig.size()){
        sscontent_json << "\"reciprocal_lattice_type_orig\":\"" << reciprocal_lattice_type_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_type_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_variation_type_orig.size()){
        sscontent_json << "\"reciprocal_lattice_variation_type_orig\":\"" << reciprocal_lattice_variation_type_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_variation_type_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_letters_orig.size()){
        vector<string> Wyckoff_letters_orig_set; 
        aurostd::string2tokens(Wyckoff_letters_orig,Wyckoff_letters_orig_set,";");
        vector<string> tmp_content;
        for(uint w=0;w<Wyckoff_letters_orig_set.size();w++){
          vector<string> Wyckoff_tokens;
          aurostd::string2tokens(Wyckoff_letters_orig_set[w],Wyckoff_tokens,",");
          tmp_content.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(Wyckoff_tokens,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_letters_orig\":[" << aurostd::joinWDelimiter(tmp_content,",") << "]" << eendl;
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_letters_orig\":\"" << Wyckoff_letters_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_letters_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_multiplicities_orig.size()){
        vector<string> Wyckoff_multiplicities_orig_set; 
        aurostd::string2tokens(Wyckoff_multiplicities_orig,Wyckoff_multiplicities_orig_set,";");
        vector<string> tmp_content;
        for(uint w=0;w<Wyckoff_multiplicities_orig_set.size();w++){
          vector<string> Wyckoff_tokens;
          aurostd::string2tokens(Wyckoff_multiplicities_orig_set[w],Wyckoff_tokens,",");
          tmp_content.push_back("["+aurostd::joinWDelimiter(Wyckoff_tokens,",")+"]");
        }
        sscontent_json << "\"Wyckoff_multiplicities_orig\":[" << aurostd::joinWDelimiter(tmp_content,",") << "]" << eendl;
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_multiplicities_orig\":\"" << Wyckoff_multiplicities_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_multiplicities_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_site_symmetries_orig.size()){
        vector<string> Wyckoff_site_symmetries_orig_set; 
        aurostd::string2tokens(Wyckoff_site_symmetries_orig,Wyckoff_site_symmetries_orig_set,";");
        vector<string> tmp_content;
        for(uint w=0;w<Wyckoff_site_symmetries_orig_set.size();w++){
          vector<string> Wyckoff_tokens;
          aurostd::string2tokens(Wyckoff_site_symmetries_orig_set[w],Wyckoff_tokens,",");
          tmp_content.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(Wyckoff_tokens,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_site_symmetries_orig\":[" << aurostd::joinWDelimiter(tmp_content,",") << "]" << eendl;
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_site_symmetries_orig\":\"" << Wyckoff_site_symmetries_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_site_symmetries_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");


      //DX20180823 - added more symmetry info - START
      // SYMMETRY
      //////////////////////////////////////////////////////////////////////////
      if(crystal_family.size()){
        sscontent_json << "\"crystal_family\":\"" << crystal_family << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_family\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(crystal_system.size()){
        sscontent_json << "\"crystal_system\":\"" << crystal_system << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_system\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(crystal_class.size()){
        sscontent_json << "\"crystal_class\":\"" << crystal_class << "\"" << eendl;
      } else {
        if(PRINT_NULL){ sscontent_json << "\"crystal_class\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Hermann_Mauguin.size()){
        sscontent_json << "\"point_group_Hermann_Mauguin\":\"" << point_group_Hermann_Mauguin << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Hermann_Mauguin\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Schoenflies.size()){
        sscontent_json << "\"point_group_Schoenflies\":\"" << point_group_Schoenflies << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Schoenflies\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_orbifold.size()){
        sscontent_json << "\"point_group_orbifold\":\"" << point_group_orbifold << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_orbifold\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_type.size()){
        sscontent_json << "\"point_group_type\":\"" << point_group_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_order!=AUROSTD_NAN){
        sscontent_json << "\"point_group_order\":" << point_group_order << eendl; //DX20190124 - changed to number, not string 
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_order\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_structure.size()){
        sscontent_json << "\"point_group_structure\":\"" << point_group_structure << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_structure\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_type.size()){
        sscontent_json << "\"Bravais_lattice_lattice_type\":\"" << Bravais_lattice_lattice_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_variation_type.size()){
        sscontent_json << "\"Bravais_lattice_lattice_variation_type\":\"" << Bravais_lattice_lattice_variation_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_variation_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_system.size()){
        sscontent_json << "\"Bravais_lattice_lattice_system\":\"" << Bravais_lattice_lattice_system << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_system\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_type.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_type\":\"" << Bravais_superlattice_lattice_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_variation_type.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_variation_type\":\"" << Bravais_superlattice_lattice_variation_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_variation_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_system.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_system\":\"" << Bravais_superlattice_lattice_system << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_system\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_superlattice.size()){
        sscontent_json << "\"Pearson_symbol_superlattice\":\"" << Pearson_symbol_superlattice << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_superlattice\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(vreciprocal_geometry.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"reciprocal_geometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vreciprocal_geometry,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_geometry\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_volume_cell!=AUROSTD_NAN) {
        sscontent_json << "\"reciprocal_volume_cell\":" << reciprocal_volume_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_volume_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_type.size()){
        sscontent_json << "\"reciprocal_lattice_type\":\"" << reciprocal_lattice_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_variation_type.size()){
        sscontent_json << "\"reciprocal_lattice_variation_type\":\"" << reciprocal_lattice_variation_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_variation_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_letters.size()){
        vector<string> Wyckoff_letters_set; 
        aurostd::string2tokens(Wyckoff_letters,Wyckoff_letters_set,";");
        vector<string> tmp_content;
        for(uint w=0;w<Wyckoff_letters_set.size();w++){
          vector<string> Wyckoff_tokens;
          aurostd::string2tokens(Wyckoff_letters_set[w],Wyckoff_tokens,",");
          tmp_content.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(Wyckoff_tokens,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_letters\":[" << aurostd::joinWDelimiter(tmp_content,",") << "]" << eendl;
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_letters\":\"" << Wyckoff_letters << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_letters\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_multiplicities.size()){
        vector<string> Wyckoff_multiplicities_set; 
        aurostd::string2tokens(Wyckoff_multiplicities,Wyckoff_multiplicities_set,";");
        vector<string> tmp_content;
        for(uint w=0;w<Wyckoff_multiplicities_set.size();w++){
          vector<string> Wyckoff_tokens;
          aurostd::string2tokens(Wyckoff_multiplicities_set[w],Wyckoff_tokens,",");
          tmp_content.push_back("["+aurostd::joinWDelimiter(Wyckoff_tokens,",")+"]");
        }
        sscontent_json << "\"Wyckoff_multiplicities\":[" << aurostd::joinWDelimiter(tmp_content,",") << "]" << eendl;
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_multiplicities\":\"" << Wyckoff_multiplicities << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_multiplicities\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_site_symmetries.size()){
        vector<string> Wyckoff_site_symmetries_set; 
        aurostd::string2tokens(Wyckoff_site_symmetries,Wyckoff_site_symmetries_set,";");
        vector<string> tmp_content;
        for(uint w=0;w<Wyckoff_site_symmetries_set.size();w++){
          vector<string> Wyckoff_tokens;
          aurostd::string2tokens(Wyckoff_site_symmetries_set[w],Wyckoff_tokens,",");
          tmp_content.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(Wyckoff_tokens,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_site_symmetries\":[" << aurostd::joinWDelimiter(tmp_content,",") << "]" << eendl;
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_site_symmetries\":\"" << Wyckoff_site_symmetries << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_site_symmetries\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //DX20190208 - added anrl info - START
      // ANRL
      //////////////////////////////////////////////////////////////////////////
      if(anrl_label_orig.size()){
        sscontent_json << "\"anrl_label_orig\":\"" << anrl_label_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"anrl_label_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(anrl_parameter_list_orig.size()){
        vector<string> anrl_parameters_vector_orig; aurostd::string2tokens(anrl_parameter_list_orig,anrl_parameters_vector_orig,",");
        sscontent_json << "\"anrl_parameter_list_orig\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(anrl_parameters_vector_orig,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"anrl_parameter_list_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(anrl_parameter_values_orig.size()){
        vector<string> anrl_values_vector_orig; aurostd::string2tokens(anrl_parameter_values_orig,anrl_values_vector_orig,",");
        sscontent_json << "\"anrl_parameter_values_orig\":[" << aurostd::joinWDelimiter(anrl_values_vector_orig,",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"anrl_parameter_values_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(anrl_label_relax.size()){
        sscontent_json << "\"anrl_label_relax\":\"" << anrl_label_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"anrl_label_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(anrl_parameter_list_relax.size()){
        vector<string> anrl_parameters_vector_relax; aurostd::string2tokens(anrl_parameter_list_relax,anrl_parameters_vector_relax,",");
        sscontent_json << "\"anrl_parameter_list_relax\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(anrl_parameters_vector_relax,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"anrl_parameter_list_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(anrl_parameter_values_relax.size()){
        vector<string> anrl_values_vector_relax; aurostd::string2tokens(anrl_parameter_values_relax,anrl_values_vector_relax,",");
        sscontent_json << "\"anrl_parameter_values_relax\":[" << aurostd::joinWDelimiter(anrl_values_vector_relax,",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"anrl_parameter_values_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //DX20190208 - added anrl info - END

      // AGL/AEL
      //////////////////////////////////////////////////////////////////////////
      if(agl_thermal_conductivity_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_thermal_conductivity_300K\":" << agl_thermal_conductivity_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_thermal_conductivity_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_debye!=AUROSTD_NAN) {
        sscontent_json << "\"agl_debye\":" << agl_debye << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_debye\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_acoustic_debye!=AUROSTD_NAN) {
        sscontent_json << "\"agl_acoustic_debye\":" << agl_acoustic_debye << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_acoustic_debye\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_gruneisen!=AUROSTD_NAN) {
        sscontent_json << "\"agl_gruneisen\":" << agl_gruneisen << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_gruneisen\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_heat_capacity_Cv_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_heat_capacity_Cv_300K\":" << agl_heat_capacity_Cv_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_heat_capacity_Cv_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_heat_capacity_Cp_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_heat_capacity_Cp_300K\":" << agl_heat_capacity_Cp_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_heat_capacity_Cp_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_thermal_expansion_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_thermal_expansion_300K\":" << agl_thermal_expansion_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_thermal_expansion_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_bulk_modulus_static_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_bulk_modulus_static_300K\":" << agl_bulk_modulus_static_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_bulk_modulus_static_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_bulk_modulus_isothermal_300K\":" << agl_bulk_modulus_isothermal_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_bulk_modulus_isothermal_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_poisson_ratio_source.size()) {
        sscontent_json << "\"agl_poisson_ratio_source\":\"" << agl_poisson_ratio_source << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_poisson_ratio_source\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_free_energy_300K_cell!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_free_energy_300K_cell\":" << agl_vibrational_free_energy_300K_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_free_energy_300K_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_free_energy_300K_atom!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_free_energy_300K_atom\":" << agl_vibrational_free_energy_300K_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_free_energy_300K_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_entropy_300K_cell!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_entropy_300K_cell\":" << agl_vibrational_entropy_300K_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_entropy_300K_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_entropy_300K_atom!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_entropy_300K_atom\":" << agl_vibrational_entropy_300K_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_entropy_300K_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_poisson_ratio!=AUROSTD_NAN) {
        sscontent_json << "\"ael_poisson_ratio\":" << ael_poisson_ratio << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_poisson_ratio\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_voigt!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_voigt\":" << ael_bulk_modulus_voigt << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_voigt\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_reuss!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_reuss\":" << ael_bulk_modulus_reuss << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_reuss\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_voigt!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_voigt\":" << ael_shear_modulus_voigt << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_voigt\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_reuss!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_reuss\":" << ael_shear_modulus_reuss << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_reuss\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_vrh\":" << ael_bulk_modulus_vrh << eendl; //CT20190117
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_vrh\":null" << eendl; //CT20190117
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_vrh\":" << ael_shear_modulus_vrh << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_vrh\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_elastic_anisotropy!=AUROSTD_NAN) { //CO20181129
        sscontent_json << "\"ael_elastic_anisotropy\":" << ael_elastic_anisotropy << eendl; //CO20181129
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_elastic_anisotropy\":null" << eendl; //CO20181129
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_youngs_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_youngs_modulus_vrh\":" << ael_youngs_modulus_vrh << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_youngs_modulus_vrh\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_speed_sound_transverse!=AUROSTD_NAN) {
        sscontent_json << "\"ael_speed_sound_transverse\":" << ael_speed_sound_transverse << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_speed_sound_transverse\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_speed_sound_longitudinal!=AUROSTD_NAN) {
        sscontent_json << "\"ael_speed_sound_longitudinal\":" << ael_speed_sound_longitudinal << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_speed_sound_longitudinal\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      ////////////////////////////////////////////////////////////////////////// 

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_speed_sound_average!=AUROSTD_NAN) {
        sscontent_json << "\"ael_speed_sound_average\":" << ael_speed_sound_average << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_speed_sound_average\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_pughs_modulus_ratio!=AUROSTD_NAN) {
        sscontent_json << "\"ael_pughs_modulus_ratio\":" << ael_pughs_modulus_ratio << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_pughs_modulus_ratio\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_debye_temperature!=AUROSTD_NAN) {
        sscontent_json << "\"ael_debye_temperature\":" << ael_debye_temperature << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_debye_temperature\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_applied_pressure!=AUROSTD_NAN) {
        sscontent_json << "\"ael_applied_pressure\":" << ael_applied_pressure << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_applied_pressure\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_average_external_pressure!=AUROSTD_NAN) {
        sscontent_json << "\"ael_average_external_pressure\":" << ael_average_external_pressure << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_average_external_pressure\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////ME20191105
      if ((ael_stiffness_tensor.rows == 6) && (ael_stiffness_tensor.cols == 6)) {
        sscontent_json << "\"ael_stiffness_tensor\":[";
        for (int i = 1; i <= 6; i++) {
          sscontent_json << "[";
          for (int j = 1; j <= 6; j++) sscontent_json << ael_stiffness_tensor[i][j] << ((j < 6)?",":"]");
          sscontent_json << ((i < 6)?",":"");
        }
        sscontent_json << "]";
      } else {
        if (PRINT_NULL) sscontent_json << "\"ael_stiffness_tensor\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////ME20191105
      if ((ael_compliance_tensor.rows == 6) && (ael_compliance_tensor.cols == 6)) {
        sscontent_json << "\"ael_compliance_tensor\":[";
        for (int i = 1; i <= 6; i++) {
          sscontent_json << "[";
          for (int j = 1; j <= 6; j++) sscontent_json << ael_compliance_tensor[i][j] << ((j < 6)?",":"]");
          sscontent_json << ((i < 6)?",":"");
        }
        sscontent_json << "]";
      } else {
        if (PRINT_NULL) sscontent_json << "\"ael_compliance_tensor\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // BADER
      //////////////////////////////////////////////////////////////////////////
      if(vbader_net_charges.size()) {
        sscontent_json << "\"bader_net_charges\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vbader_net_charges,6),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"bader_net_charges\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vbader_atomic_volumes.size()) {
        sscontent_json << "\"bader_atomic_volumes\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vbader_atomic_volumes,4),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"bader_atomic_volumes\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // FILES
      //////////////////////////////////////////////////////////////////////////
      if(vfiles.size()) {
        aurostd::sort(vfiles);
        sscontent_json << "\"files\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vfiles,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"files\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // CPUS
      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_Model.size()) {
        sscontent_json << "\"node_CPU_Model\":\"" << node_CPU_Model << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_Model\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_Cores!=AUROSTD_NAN) {
        sscontent_json << "\"node_CPU_Cores\":" << node_CPU_Cores << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_Cores\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_MHz!=AUROSTD_NAN) {
        sscontent_json << "\"node_CPU_MHz\":" << node_CPU_MHz << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_MHz\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_RAM_GB!=INF) {
        sscontent_json << "\"node_RAM_GB\":" << node_RAM_GB << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_RAM_GB\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // VERSION/DATE
      //////////////////////////////////////////////////////////////////////////
      if(aflow_version.size()) {
        sscontent_json << "\"aflow_version\":\"" << aflow_version << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_version\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(catalog.size()) {
        sscontent_json << "\"catalog\":\"" << catalog << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"catalog\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      sscontent_json << "\"aflowlib_version\":\"" << string(AFLOW_VERSION) << "\"" << eendl;  //CO20170613
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //[CO20200624 - OBSOLETE]sscontent_json << "\"aflowlib_date\":\"" << aurostd::get_datetime() << "_GMT-5\"" << eendl;
      //[CO20200624 - OBSOLETE]vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      if(vaflowlib_date.size()) {
        sscontent_json << "\"aflowlib_date\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vaflowlib_date,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflowlib_date\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
      vcontent_json.clear();

      sss << endl;
    } // json

    return sss.str();
  }

  //  bool aflowlib2file(
  string _aflowlib_entry::aflowlib2file(string file,string mode) {
    string aflowlib_out=aflowlib2string(mode);
    aurostd::string2file(aflowlib_out,file);
    return aflowlib_out;
  }
}

//CO20171202 - apennsy fixes!
namespace aflowlib {
  void _aflowlib_entry::correctBadDatabase(bool verbose,ostream& oss){
    ofstream FileMESSAGE;
    return correctBadDatabase(FileMESSAGE,verbose,oss);
  }

  void _aflowlib_entry::correctBadDatabase(ofstream& FileMESSAGE,bool verbose,ostream& oss){
    //CO20180828 - LIB2 also contains unaries //so far we only know of bad binaries
    //APENNSY neglect - LIB2 only //CO20180828 - LIB2 also contains unaries  //binaries only
    string soliloquy = XPID + "_aflowlib_entry::correctBadDatabase():";
    stringstream message;
    if(vspecies_pp.size()==1 || vspecies_pp.size()==2) {
      string pseudoA="",pseudoB="";
      pseudoA = vspecies_pp[0];
      if(vspecies_pp.size()==2){pseudoB = vspecies_pp[1];}
      //[OBSOLETE CO20180828]string pseudoA = vspecies_pp[0];
      //[OBSOLETE CO20180828]string pseudoB = vspecies_pp[1];
      // tiny corrections
      //gamma_IrV
      if(pseudoA=="Cd" && pseudoB=="Pt" && prototype=="181") {
        enthalpy_formation_atom -= 0.0013;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        if(verbose){
          message << "Fixing enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      //gamma_IrV
      if(pseudoA=="Ir" && pseudoB=="V_sv" && prototype=="291") {
        enthalpy_formation_cell -= 0.001;
        enthalpy_formation_atom -= 0.0005;
        enthalpy_cell -= 0.001;
        enthalpy_atom -= 0.005;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // HfPd
      if(pseudoA=="Hf_pv" && pseudoB=="Pd_pv" && prototype=="192") {
        enthalpy_formation_atom -= 0.003;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        enthalpy_atom = enthalpy_formation_atom;
        enthalpy_cell = natoms * enthalpy_atom;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // sigma
      if(pseudoA=="Ir" && pseudoB=="Nb_sv" && prototype=="600.ABBAB") {
        enthalpy_formation_cell += 0.001;
        enthalpy_formation_atom += 0.0005;
        enthalpy_cell += 0.001;
        enthalpy_atom += 0.005;
        enthalpy_formation_cell += 0.001;
        enthalpy_formation_atom += 0.0005;
        enthalpy_cell += 0.001;
        enthalpy_atom += 0.005;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // sigma
      if(pseudoA=="Os_pv" && pseudoB=="Re_pv" && prototype=="122") {
        enthalpy_formation_atom -= 0.001;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        enthalpy_atom = enthalpy_formation_atom;
        enthalpy_cell = natoms * enthalpy_atom;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
    }
  }
  bool _aflowlib_entry::ignoreBadDatabase() const{
    string reason;
    return ignoreBadDatabase(reason);
  }
  bool _aflowlib_entry::ignoreBadDatabase(string& reason) const{
    reason="";
    //grab bad database protos
    vector<string> _DEVIL_PROTOTYPES_;
    aurostd::string2tokens(_DEVIL_PROTOTYPES_STRING_,_DEVIL_PROTOTYPES_,",");

    //so far we only know of bad binaries
    //we need something more robust than just exact string match, case: 549 and 549.bis vs. 549.tetra
    bool match=false;
    //DEVIL
    for(uint di=0;di<_DEVIL_PROTOTYPES_.size() && !match;di++){if(pflow::prototypeMatch(prototype, _DEVIL_PROTOTYPES_[di])){match=true;}}
    if(match){
      reason=compound+":"+prototype+" is ill-calculated in the database";
      return true;
    }
    //find .old's
    if(1){
      uint prototype_size=prototype.size();
      string search_string=".old";uint search_string_size=search_string.size();
      if(prototype_size>search_string_size && prototype.substr(prototype_size-search_string_size,search_string_size)==search_string){  //look only at the end of the prototype
        reason=compound+":"+prototype+" is ill-calculated in the database";
        return true;
      }
    }
    //APENNSY neglect - LIB2 only //CO20180828 - LIB2 also contains unaries  //binaries only
    if(vspecies_pp.size()==1 || vspecies_pp.size()==2) {
      string pseudoA="",pseudoB="";
      pseudoA = vspecies_pp[0];
      if(vspecies_pp.size()==2){pseudoB = vspecies_pp[1];}
      //[OBSOLETE CO20180828]string pseudoA = vspecies_pp[0];
      //[OBSOLETE CO20180828]string pseudoB = vspecies_pp[1];
      // bad Ag is a wrong relaxation
      if((pseudoA=="Ag" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ag" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ag is a wrong relaxation
      if((pseudoA=="Ag" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Ag" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Au is a wrong relaxation
      if((pseudoA=="Au" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Au" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Al_h" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Al_h" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Al_h" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ca_sv is a wrong relaxation
      if((pseudoA=="Ca_sv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ca_sv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ca_sv is a wrong relaxation
      if((pseudoA=="Ca_sv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Ca_sv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cd is a wrong relaxation
      if((pseudoA=="Cd" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Cd" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cu_pv is a wrong relaxation
      if((pseudoA=="Cu_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Cu_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cu_pv is a wrong relaxation
      if((pseudoA=="Cu_pv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Cu_pv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Fe_pv is a wrong relaxation
      if((pseudoA=="Fe_pv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Fe_pv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Fe_pv is a wrong relaxation
      if((pseudoA=="Fe_pv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Fe_pv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ge_h is a wrong relaxation
      if((pseudoA=="Ge_h" && pflow::prototypeMatch(prototype,"305")) || (pseudoB=="Ge_h" && pflow::prototypeMatch(prototype,"306"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad In_d is a wrong relaxation
      if((pseudoA=="In_d" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="In_d" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ir is a wrong relaxation
      if((pseudoA=="Ir" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ir" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad K_sv is a wrong relaxation
      if((pseudoA=="K_sv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="K_sv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad K_sv is a wrong relaxation
      if((pseudoA=="K_sv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="K_sv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad La is a wrong relaxation
      if((pseudoA=="La" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="La" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad La is a wrong relaxation
      if((pseudoA=="La" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="La" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Li_sv is a wrong relaxation
      if((pseudoA=="Li_sv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Li_sv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Li_sv is a wrong relaxation
      if((pseudoA=="Li_sv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Li_sv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Na_pv is a wrong relaxation
      if((pseudoA=="Na_pv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Na_pv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Na_pv is a wrong relaxation
      if((pseudoA=="Na_pv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Na_pv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ni_pv is a wrong relaxation
      if((pseudoA=="Ni_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ni_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ni_pv is a wrong relaxation
      if((pseudoA=="Ni_pv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Ni_pv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pb_d is a wrong relaxation
      if((pseudoA=="Pb_d" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Pb_d" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pb_d is a wrong relaxation
      if((pseudoA=="Pb_d" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Pb_d" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pd_pv is a wrong relaxation
      if((pseudoA=="Pd_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Pd_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pd_pv is a wrong relaxation
      if((pseudoA=="Pd_pv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Pd_pv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pt is a wrong relaxation
      if((pseudoA=="Pt" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Pt" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pt is a wrong relaxation
      if((pseudoA=="Pt" && pflow::prototypeMatch(prototype,"317")) || (pseudoB=="Pt" && pflow::prototypeMatch(prototype,"318"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Rh_pv is a wrong relaxation
      if((pseudoA=="Rh_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Rh_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"305")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"306"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ta_pv is a wrong relaxation
      if((pseudoA=="Ta_pv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Ta_pv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ta_pv is a wrong relaxation
      if((pseudoA=="Ta_pv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Ta_pv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad B_h is a wrong relaxation
      if((pseudoA=="B_h" && pflow::prototypeMatch(prototype,"317")) || (pseudoB=="B_h" && pflow::prototypeMatch(prototype,"318"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }

      // sigma
      if(pseudoA=="Os_pv" && pseudoB=="Re_pv" && pflow::prototypeMatch(prototype,"448")) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // wrong channel, bug
      if(pseudoA=="Rh_pv" && pseudoB=="Zr_sv" && pflow::prototypeMatch(prototype,"381")) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
    }
    return false;
  }
} // namespace aflowlib

namespace aflowlib {
  string _aflowlib_entry::getPathAURL(ostream& oss, bool load_from_common){ //CO20200404
    ofstream FileMESSAGE;
    return getPathAURL(FileMESSAGE, oss, load_from_common);
  }
  string _aflowlib_entry::getPathAURL(ofstream& FileMESSAGE,ostream& oss, bool load_from_common){  //CO20200404
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "_aflowlib_entry::getPathAURL():";
    stringstream message;
    string path = "";
    if (aurl.empty()) {return path;}
    vector<string> tokens;
    aurostd::string2tokens(aurl, tokens, ":");
    //LIB1 presents problems here (3 colons): aflowlib.duke.edu:AFLOWDATA/LIB1_RAW/Pt:PAW_PBE:05Jan2001/A6
    if(0){
      if (tokens.size() != 2) {
        message << "Odd AURL format for entry " << auid << ": " << aurl;
        pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, FileMESSAGE, oss, _LOGGER_WARNING_);
        return path;
      }
    }

    //instead just erase first item, join others, assume we're okay...
    tokens.erase(tokens.begin());
    path=aurostd::joinWDelimiter(tokens,":");

    string server="",path_full="";
    if(1&&load_from_common){
      //attempt 1: try replacing _RAW with _LIB
      if(1){
        server="/www";
        path_full=path;
        aurostd::StringSubst(path_full,"_RAW","_LIB");
        path_full=server+"/"+path_full;
        if(LDEBUG){cerr << soliloquy << " attempt 1 path=" << path_full << endl;}
        if(aurostd::IsDirectory(path_full)){return path_full;}
      }

      //attempt 2: try finding LIB directory
      if(1){
        server="/common";
        path_full=path;
        aurostd::StringSubst(path_full,"AFLOWDATA/","");
        aurostd::StringSubst(path_full,"ICSD_WEB","ICSD/LIB"); //CO20200223
        aurostd::StringSubst(path_full,"_RAW","/LIB");
        path_full=server+"/"+path_full;
        if(LDEBUG){cerr << soliloquy << " attempt 2 path=" << path_full << endl;}
        if(aurostd::IsDirectory(path_full)){return path_full;}
      }

      //attempt 3: try no replacement (RAW)
      if(1){
        server="/www";
        path_full=server+"/"+path;
        if(LDEBUG){cerr << soliloquy << " attempt 3 path=" << path_full << endl;}
        if(aurostd::IsDirectory(path_full)){return path_full;}
      }
    }
    if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")){server=XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");}
    else{server="aflowlib.duke.edu";}

    path_full=server+"/"+path;
    return path_full;
  }
}

// **************************************************************************
// directory2auid
// auid2directory
// auid2present
// **************************************************************************
namespace aflowlib {
  bool _aflowlib_entry::directory2auid(string directory) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "_aflowlib_entry::directory2auid: BEGIN" << endl;
    auid="";
    vauid.clear();

    bool conflict=TRUE; 
    while (conflict) {
      uint64_t crc=0;
      // DONT TOUCH THE AUID FLOW

      vector<string> vfiles2;
      for(uint iext=1;iext<XHOST.vext.size();iext++) {  // SKIP uncompressed
        vfiles2.push_back("OUTCAR.relax1"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.relax2"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.relax3"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.relax4"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.static"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.bands"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR"+XHOST.vext.at(iext));
      }

      if(LDEBUG) cerr << XPID << "_aflowlib_entry::directory2auid: [0]" << endl;

      for(uint i=0;i<vfiles2.size();i++) 
        if(aurostd::FileExist(directory+"/"+vfiles2.at(i)))
          crc=aurostd::crc64(crc,aurostd::efile2string(directory+"/"+vfiles2.at(i))); // DONT TOUCH THIS
      auid="aflow:"+aurostd::crc2string(crc);

      if(LDEBUG) cerr << XPID << "_aflowlib_entry::directory2auid: [1]" << endl;

      if(LDEBUG) cerr << XPID << "_aflowlib_entry::directory2auid: auid=" << auid << endl;
      conflict=FALSE;
      string aurl_found;
      if(aflowlib::auid2present(auid,aurl_found,1)) {
        if(LDEBUG) cerr << XPID << "_aflowlib_entry::directory2auid: conflict auid=" << auid << endl;	
        cerr << XPID << "[WARNING]  _aflowlib_entry::directory2auid: CONFLICT POTENTIAL " << auid << " " << aurl_found << " " << aurl << endl;
        if(aurl_found!=aurl) { // avoid conflict with yourself
          string salt="AUID_salt["+aurostd::utype2string<long double>(aurostd::get_useconds())+"]";
          cerr << XPID << "[WARNING]  _aflowlib_entry::directory2auid: CONFLICT TRUE      " << auid << " " << aurl_found << " " << aurl << "  " << salt << endl;
          string file=vfiles2.at(0);
          //
          for(uint iext=0;iext<XHOST.vext.size();iext++) { aurostd::StringSubst(file,XHOST.vext.at(iext),""); }
          stringstream sss;aurostd::efile2stringstream(directory+"/"+file+DEFAULT_KZIP_EXT,sss); sss << endl << salt << endl;
          aurostd::execute("mv \""+directory+"/"+file+DEFAULT_KZIP_EXT+"\""+" \""+directory+"/"+file+".conflict_auid"+DEFAULT_KZIP_EXT+"\"");
          aurostd::stringstream2compressfile(DEFAULT_KZIP_BIN,sss,directory+"/"+file);
          //
          conflict=TRUE; // recheck
        } else {
          cerr << XPID << "[WARNING]  _aflowlib_entry::directory2auid: CONFLICT TRIVIAL   " << auid << " " << aurl_found << " " << aurl << endl;
        }
      }
    }

    // now it has auid   
    vauid.clear();
    aflowlib::auid2vauid(auid,vauid);

    if(LDEBUG) cerr << "directory2auid: END" << endl;

    cout << XPID << "_aflowlib_entry::directory2auid: DIRECTORY=" << directory << endl; // DONT TOUCH THIS
    cout << XPID << "_aflowlib_entry::directory2auid: AURL_ID=" << aurostd::crc2string(aurostd::crc64(0,directory)) << endl; // DONT TOUCH THIS

    return TRUE;
  }

  bool json2aflowlib(const string& json,string key,string& value) { //SC20200415
    // return TRUE if something has been found
    value="";
    key="\""+key+"\":";
    string::size_type start,end;
    start=json.find(key);
    if(start!=string::npos) {
      start+=key.length();
      end=json.find("\":",start);
      if(end!=string::npos){
        value=json.substr(start,end-start);
        end=value.find_last_of(",");
        value=value.substr(0,end);
      } else {
        end=json.find("}",start);
        value=json.substr(start,end-start);
      }
      //    if((value[0]=='\"') && (value[value.size()-1]=='\"')) value=value.substr(1,value.size()-2);  // Remove quotes
    } else {
      value="";
    }
    // cleanup
    aurostd::StringSubst(value,"[","");  // Remove brakets
    aurostd::StringSubst(value,"]","");  // Remove brakets
    aurostd::StringSubst(value,"\"",""); // Remove quotes

    return !value.empty();
  }

  uint auid2present(string auid,string& aurl,int mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::auid2present: BEGIN mode=" << mode << endl;
    string loop="",json="";aurl="";
    if(auid=="" || auid.size()!=22) { cerr << XPID << "aflowlib::auid2present: auid.size() needs to be 22 characters long" << endl; return FALSE;}

    // [OBSOLETE]    if(mode==0) {
    // [OBSOLETE]    if(XHOST_vAUID.size()==0 || XHOST_vAURL.size()==0) init::InitGlobalObject("vLIBS","",FALSE);
    // [OBSOLETE]    if(LDEBUG) cerr << XPID << "aflowlib::auid2present: [4] XHOST_vAURL.size()=" << XHOST_vAURL.size() << "  XHOST_vAUID.size()=" << XHOST_vAUID.size() << endl;
    // [OBSOLETE]    bool found=FALSE;
    // [OBSOLETE]    for(uint j=0;j<XHOST_vAUID.size()&&!found;j++) {
    // [OBSOLETE]  	if(LDEBUG && XHOST_vAUID.at(j)==auid) cerr << "[" << auid << "] [" << XHOST_vAUID.at(j) << "] [" << XHOST_vAURL.at(j) << "]" << " [" << j << "]" << endl;
    // [OBSOLETE]  	if(XHOST_vAUID.at(j)==auid) {found=TRUE;aurl=XHOST_vAURL.at(j);}
    // [OBSOLETE]    }
    // [OBSOLETE]    if(LDEBUG) cerr << XPID << "aflowlib::auid2present: END  auid=" << auid << "  aurl=" << aurl << endl;
    // [OBSOLETE]    return found;
    // [OBSOLETE]   }
    if(mode==1) { // PICK THIS ONE DEFAULT
      //  bool aflowlib::auid2present(string auid,string& aurl)
      string jsonl_file=XHOST_LIBRARY_JSONL+"/"; for(uint i=0;i<8;i++) jsonl_file+=auid.at(i); jsonl_file+=".jsonl";
      bool found=FALSE;
      jsonl_file=aurostd::CleanFileName(jsonl_file);
      for(uint i=0;i<XHOST.vext.size()&&aurl.empty()&&!found;i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::auid2present: TESTING=" << jsonl_file << XHOST.vext.at(i) << endl; 
        //	cout << XPID << "aflowlib::auid2present: TESTING=" << jsonl_file << XHOST.vext.at(i) << endl; 
        if(aurostd::FileExist(jsonl_file+XHOST.vext.at(i))) {
          found=TRUE;
          if(LDEBUG) cerr << XPID << "aflowlib::auid2present: FOUND=" << jsonl_file << XHOST.vext.at(i) << endl; 
          //  cout << XPID << "aflowlib::auid2present: FOUND=" << jsonl_file << XHOST.vext.at(i) << endl; 
          json=aurostd::execute2string(XHOST.vcat.at(i)+" "+jsonl_file+XHOST.vext.at(i)+" | grep "+auid);
          aflowlib::json2aflowlib(json,"aurl",aurl);
          aflowlib::json2aflowlib(json,"loop",loop);
        }
      }
      if(LDEBUG) cerr << XPID << "aflowlib::auid2present: END  auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      cout << XPID << "aflowlib::auid2present: auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      return json.size();
    }
    if(mode==2) { // not that faster and does not keep an outside vAUID table so it does not see the TRIVIAL CONFLICTS
      //  bool aflowlib::auid2present(string auid,string& aurl)
      string jsonl_file=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_AUID)+"/"+aflowlib::auid2directory(auid)+"/RAW/aflowlib.json";
      bool found=FALSE;
      jsonl_file=aurostd::CleanFileName(jsonl_file);
      for(uint i=0;i<XHOST.vext.size()&&!found;i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::auid2present: TESTING=" << jsonl_file << XHOST.vext.at(i) << endl; 
        if(aurostd::FileExist(jsonl_file+XHOST.vext.at(i))) {
          found=TRUE;
          if(LDEBUG) cerr << XPID << "aflowlib::auid2present: FOUND=" << jsonl_file << XHOST.vext.at(i) << endl;
          json=aurostd::execute2string(XHOST.vcat.at(i)+" "+jsonl_file+XHOST.vext.at(i));
          aflowlib::json2aflowlib(json,"aurl",aurl);
          aflowlib::json2aflowlib(json,"loop",loop);
        }
      }
      if(LDEBUG) cerr << XPID << "aflowlib::auid2present: END  auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      cout << XPID << "aflowlib::auid2present: auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      return json.size();
    }

    return FALSE;
  }

  // [OBSOLETE]   time ./aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/Cu_pvHgSn/TFCC004.CAB
  // [OBSOLETE]   time aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/Cu_pvHgSn/TFCC004.CAB  
  // [OBSOLETE]   time ./aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/AgCdCo/TFCC001.ABC   
  // [OBSOLETE]   time aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/AgCdCo/TFCC001.ABC  
  // [OBSOLETE]   time ./aflow --beep --force --showPID --lib2raw="/common/LIB4/LIB/AgCdCoZr_sv:PAW_PBE/ABCD_cF16_216_c_d_b_a.ABCD"
  // [OBSOLETE]   time aflow --beep --force --showPID --lib2raw="/common/LIB4/LIB/AgCdCoZr_sv:PAW_PBE/ABCD_cF16_216_c_d_b_a.ABCD"


  uint auid2vauid(const string auid, deque<string>& vauid) {                // splits the auid into vauid
    vauid.clear();
    //    vauid.push_back(auid.substr(0,6)); for(uint i=6;i<=20;i+=2) vauid.push_back(auid.substr(i,2));  // splitting aflow:/ab/cd..
    vauid.push_back(auid.substr(0,8)); for(uint i=8;i<=20;i+=2) vauid.push_back(auid.substr(i,2));  // splitting aflow:ab/cd..
    return vauid.size();
  }

  string auid2directory(const string auid) {                                // gives AUID directory from existence of vauid
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::auid2directory: BEGIN" << endl;
    string directory;
    deque<string> vauid;
    aflowlib::auid2vauid(auid,vauid);
    if(vauid.size()>0) {
      directory=vauid.at(0);
      for(uint i=1;i<vauid.size();i++) {
        directory+="/"+vauid.at(i);
      }
    }
    if(LDEBUG) cerr << XPID << "aflowlib::auid2directory: END" << endl;
    return directory;
  }

}

// **************************************************************************
// Operate on CIFs for JMOL for one structure
// **************************************************************************
namespace aflowlib {
  bool cif2data(string file,double& a,double& b,double& c,double& alpha,double& beta,double& gamma) {
    vector<string> vline,tokens;
    aurostd::file2vectorstring(file,vline);
    for(uint i=0;i<vline.size();i++) {
      aurostd::string2tokens(vline.at(i),tokens," ");
      if(aurostd::substring2bool(vline.at(i),"_cell_length_a") && tokens.size()>1) a=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_length_b") && tokens.size()>1) b=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_length_c") && tokens.size()>1) c=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_alpha") && tokens.size()>1) alpha=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_beta") && tokens.size()>1) beta=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_gamma") && tokens.size()>1) gamma=aurostd::string2utype<double>(tokens.at(1));
    }
    return TRUE;
  }

  bool cif2oss(string file,string label,string sgnumber,ostream& oss) {
    double aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF;
    aflowlib::cif2data(file,aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
    oss << " font echo 14;color echo white; ";
    oss << " set echo 3% 97%; echo \\\"AFLOW-JSmol consortium (AFLOW V" << string(AFLOW_VERSION) << ") | entry="<< label << "  |  | ";
    if(aurostd::string2utype<int>(sgnumber)>0) {
      oss <<"  Spacegroup = "<< GetSpaceGroupName(aurostd::string2utype<int>(sgnumber)) << " (#" <<  sgnumber << ")   |";
    }
    oss << " a=" << aCIF <<"\u212B, b=" << bCIF <<"\u212B, c=" << cCIF <<"\u212B    |" << " \u03B1=" << alphaCIF <<"\u00B0, \u03B2=" << betaCIF <<"\u00B0, \u03B3=" << gammaCIF <<"\u00B0   \\\"; ";

    return TRUE;
  }
}

//BEGIN JJPR
// **************************************************************************
// GET BADER jvxl file
// **************************************************************************

namespace aflowlib {
  bool iso2oss(string file,string label, string element,string cutoff,int index,ostream& oss) {
    vector<string> kk;
    kk.resize(9);
    kk[0]="red";
    kk[1]="green";
    kk[2]="yellow";
    kk[3]="blue";
    kk[4]="orange";
    kk[5]="white";
    kk[6]="purple";
    kk[7]="brown";
    kk[8]="pink";
    oss <<"  ISOSURFACE  " << element << " \\\""<< file+"/"+label+"_Bader_"+cutoff+"_"+element+".jvxl\\\"  ;isosurface mesh;  color isosurface " << kk[index] << " translucent;";
    oss << "x = load(\\\""<< file+"/"+label+"_abader.out" << "\\\"); charges = x.split(\\\"=\\\")[2].split(\\\"(\\\")[1].split(\\\",\\\"); {*}.label = charges; label \%[label];";
    // FOR LOCAL TEST: oss <<"  ISOSURFACE  " << element << " "<< "Bader_"+cutoff+"_"+element+".jvxl  ;isosurface mesh;  color isosurface " << kk[index] << " translucent;";
    return TRUE;
  }
}
//END JJPR

// **************************************************************************
// GET SPACE GROUP for one structure
// **************************************************************************
namespace aflowlib {
  uint SGtoNSG(string sgroup) {
    string::size_type idx1;
    string strsub("#");
    idx1=sgroup.find(strsub);
    if(idx1!=string::npos)  
      return (int) atoi(sgroup.substr(sgroup.find(strsub)+strsub.length()).c_str());
    else return 0;
  }

  void _aflowlib_entry::GetSGROUP(string aflowlibentry) {
    vector<string> vaflowlib_entry;
    aurostd::string2tokens(aflowlibentry,vaflowlib_entry,"|");

    bool VERBOSE_LOCAL=(FALSE || XHOST.DEBUG);
    if(vaflowlib_entry.size()==0) {cerr << "ERROR - aflowlib_entry::GetSGROUP(): " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " file not found " << endl;exit(0);}
    vsgroup.clear();vNsgroup.clear();

    string data_aurl="";vector<string> tokens;
    // CHECK FOR HOLES
    if(XGNDSTATE_HOLES==0)  // SAFE NO HOLES IN THE XMATRIX
      if(vaflowlib_entry.size()==0) {cerr << "ERROR - aflowlib_entry::GetSGROUP(): " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " file not found " << endl;exit(0);}
    if(XGNDSTATE_HOLES==1)  // ALLOW HOLES WITH FAKE VALUES
      if(vaflowlib_entry.size()==0) {
        //cerr << "FOUND SG HOLE = " << alloy_dir << "/" << params.structures[2].at(str_number).name << "   === " << str_number << " " << structures_number[str_number]<< endl;
        vsgroup.clear();vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);
        vNsgroup.clear();vNsgroup.push_back(0);vNsgroup.push_back(0);vNsgroup.push_back(0);
        return;
      }

    for(uint i=0;i<vaflowlib_entry.size();i++) {
      if(aurostd::substring2bool(vaflowlib_entry.at(i),"aurl=")) data_aurl=aurostd::substring2string(vaflowlib_entry.at(i),"arl=");
      if(aurostd::substring2bool(vaflowlib_entry.at(i),"sg=")) {
        aurostd::string2tokens(vaflowlib_entry.at(i),tokens,",");
        if(tokens.size()==0) {cerr << "ERROR - aflowlib_entry::GetSGROUP(): geometry not enough tokens " << endl;exit(0);}
        if(tokens.size()==3) { // ok
          vsgroup.clear();
          aurostd::StringSubst(tokens.at(0),"sg=","");aurostd::StringSubst(tokens.at(0)," ","");aurostd::StringSubst(tokens.at(0),"#"," #");vsgroup.push_back(tokens.at(0));
          aurostd::StringSubst(tokens.at(1)," ","");aurostd::StringSubst(tokens.at(1),"#"," #");vsgroup.push_back(tokens.at(1));
          aurostd::StringSubst(tokens.at(2)," ","");aurostd::StringSubst(tokens.at(2),"#"," #");vsgroup.push_back(tokens.at(2));
          if(VERBOSE_LOCAL) cout << "[5] "<<"["<<tokens.at(0) << "," << tokens.at(1) << "," << tokens.at(2) <<"]" << endl;
        } else {
          vsgroup.clear();vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);
          vNsgroup.clear();vNsgroup.push_back(0);vNsgroup.push_back(0);vNsgroup.push_back(0);
        }
      }
    }
    // done, now makes the numbers
    for(uint ii=0;ii<vsgroup.size();ii++) vNsgroup.push_back(SGtoNSG(vsgroup.at(ii)));
    // if(NsgroupPRE==0) {sgroupPRE=sgroupMID;NsgroupPRE=SGtoNSG(sgroupPRE);}	
    // if(NsgroupPRE==0) {sgroupPRE=sgroupPOST;NsgroupPRE=SGtoNSG(sgroupPRE);}	
    for(uint i=vsgroup.size()-2;i<vsgroup.size();i++) {
      if(vNsgroup.at(i)==0) {
        vsgroup.at(i)=vsgroup.at(i-1);
        vNsgroup.at(i)=SGtoNSG(vsgroup.at(i));
      }
    }

    if(vNsgroup.size()!=3) for(uint ii=vNsgroup.size();ii<3;ii++) vNsgroup.push_back(0);
    if(vsgroup.size()!=3) for(uint ii=vsgroup.size();ii<3;ii++) vsgroup.push_back("NNN #0");

    return;
  }
}


// ***************************************************************************
// aflowlib::AflowlibLocator
// ***************************************************************************
namespace aflowlib { // move to web interface
  bool AflowlibLocator(const string& in,string& out,const string& mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: BEGIN" << endl;

    if(mode!="AFLOWLIB_AUID2AURL" && mode!="AFLOWLIB_AURL2AUID" && mode!="AFLOWLIB_AUID2LOOP" && mode!="AFLOWLIB_AURL2LOOP") {
      cerr << "ERROR - aflowlib::AflowlibLocator: wrong mode=" << mode << endl;
      exit(0);
    }
    if(XHOST_vAUID.size()==0) { init::InitGlobalObject("vLIBS"); }
    if(XHOST_vAURL.size()==0) { init::InitGlobalObject("vLIBS"); }
    if(XHOST_vLOOP.size()==0) { init::InitGlobalObject("vLIBS"); }
    if(XHOST_vAUID.size()!=XHOST_vAURL.size() || XHOST_vAUID.size()!=XHOST_vLOOP.size()) {
      cerr << "ERROR - aflowlib::AflowlibLocator: XHOST_vAUID.size()!=XHOST_vAURL.size() || XHOST_vAUID.size()!=XHOST_vLOOP.size()" << endl;
      cerr << "                                   XHOST_vAUID.size()=" << XHOST_vAUID.size() << endl;
      cerr << "                                   XHOST_vAURL.size()=" << XHOST_vAURL.size() << endl;
      cerr << "                                   XHOST_vLOOP.size()=" << XHOST_vLOOP.size() << endl;
      exit(0);
    }
    out="";
    for(uint i=0;i<XHOST_vAUID.size()&&out.empty();i++) {
      if(mode=="AFLOWLIB_AUID2AURL" && XHOST_vAUID.at(i)==in) out=XHOST_vAURL.at(i);
      if(mode=="AFLOWLIB_AURL2AUID" && XHOST_vAURL.at(i)==in) out=XHOST_vAUID.at(i);
      if(mode=="AFLOWLIB_AUID2LOOP" && XHOST_vAUID.at(i)==in) out=XHOST_vLOOP.at(i);
      if(mode=="AFLOWLIB_AURL2LOOP" && XHOST_vAURL.at(i)==in) out=XHOST_vLOOP.at(i);
      //     cerr << i << endl;
    }
    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: END" << endl;
    return !out.empty();
  }
} // namespace aflowlib

namespace aflowlib {
  string AflowlibLocator(string options, string mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: BEGIN" << endl;
    if(mode!="AFLOWLIB_AUID2AURL" && mode!="AFLOWLIB_AURL2AUID" && mode!="AFLOWLIB_AUID2LOOP" && mode!="AFLOWLIB_AURL2LOOP") {
      cerr << "error - aflowlib::AflowlibLocator: wrong mode=" << mode << endl;
      exit(0);
    }
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==0) {
      if(mode=="AFLOWLIB_AUID2AURL") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_auid2aurl=auid1,auid2....");
      if(mode=="AFLOWLIB_AURL2AUID") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_aurl2auid=aurl1,aurl2....");
      if(mode=="AFLOWLIB_AUID2LOOP") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_auid2loop=auid1,auid2....");
      if(mode=="AFLOWLIB_AURL2LOOP") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_aurl2loop=aurl1,aurl2....");
      exit(0);
    } 
    // move on
    stringstream output;
    string locator;
    for(uint i=0;i<tokens.size();i++) {
      if(aflowlib::AflowlibLocator(tokens.at(i),locator,mode)) {
        output << locator << endl;
      } else {
        output << tokens.at(i) << " not found" << endl;
      }
    }

    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: END" << endl;
    return output.str();
  }
} // namespace aflowlib

//AFLUX integration
//FR+CO20180329
namespace aflowlib {
  bool APIget::establish(){
    struct hostent * host = gethostbyname( Domain.c_str() );

    //[CO20181226 - OBSOLETE]PORT=80;  //CO20180401

    if ( (host == NULL) || (host->h_addr == NULL) ) {
      cerr << "Error retrieving DNS information." << endl;
      return false;
      //exit(1);
    }

    bzero(&client, sizeof(client));
    client.sin_family = AF_INET;
    client.sin_port = htons( PORT );
    memcpy(&client.sin_addr, host->h_addr, host->h_length);

    sock = socket(AF_INET, SOCK_STREAM, 0);

    if (sock < 0) {
      cerr << "Error creating socket." << endl;
      return false;
      //exit(1);
    }

    if ( connect(sock, (struct sockaddr *)&client, sizeof(client)) < 0 ) {
      close(sock);
      cerr << "Could not connect" << endl;
      return false;
      //exit(1);
    }

    stringstream ss;
    ss << "GET " << API_Path << Summons << " HTTP/1.0\r\n" ;
    //    cerr << "GET " << API_Path << Summons << " HTTP/1.0\r\n" ;
    ss << "HOST: " << Domain << "\r\n";
    ss << "Connection: close\r\n";
    ss << "\r\n";
    string request = ss.str();

    if (send(sock, request.c_str(), request.length(), 0) != (int)request.length()) {
      cerr << "Error sending request." << endl;
      return false;
      //exit(1);
    }
    return true;
  }
  void APIget::reset( string a_Summons, string a_API_Path, string a_Domain ) {
    if( a_Summons == "#" ) {
      Summons = "";
      API_Path = "/search/API/?";
      Domain = "aflowlib.duke.edu";
    } else {
      Summons = a_Summons;
      if( ! a_API_Path.empty() ) API_Path = a_API_Path;
      if( ! a_Domain.empty() ) Domain = a_Domain;
    }
  }
  ostream& operator<<( ostream& output, APIget& a ) { 
    char cur;
    bool responsedata = false;
    bool waslinefeed = false;
    if( a.establish() ) {
      while ( ! responsedata ) { //discard headers
        read(a.sock, &cur, 1);
        //cerr << cur << ":" << (int)cur << endl;
        if( waslinefeed  && cur == '\r') responsedata = true;
        if( cur == '\n' ) waslinefeed = true;
        else waslinefeed = false;
      };
      read(a.sock, &cur, 1); //discard final \n in header \r\n\r\n
      while ( read(a.sock, &cur, 1) > 0 ) output << cur; //cout << cur;
      close(a.sock);
    }
    return output;
  }
}

//DX+FR20190206 - AFLUX functionality via command line - START
// ***************************************************************************
namespace aflowlib {
  string AFLUXCall(aurostd::xoption& vpflow){

    // Performs AFLUX call based on summons input from command line

    string usage="aflow --aflux=<summons>";
    string options="";

    if(vpflow.flag("AFLUX::USAGE")) {
      //[CO20200624 - OBSOLETE]stringstream ss_usage;
      //[CO20200624 - OBSOLETE]init::ErrorOption(ss_usage,vpflow.getattachedscheme("AFLUX"),"aflowlib::AFLUXCall()",aurostd::liststring2string(usage,options));
      //[CO20200624 - OBSOLETE]return ss_usage.str();
      init::ErrorOption(vpflow.getattachedscheme("AFLUX"),"aflowlib::AFLUXCall()",aurostd::liststring2string(usage,options));
      return "";
    }

    string summons = "";
    if(vpflow.flag("AFLUX")) {
      summons=vpflow.getattachedscheme("AFLUX");
      // check if string is enclosed in double or single quotes 
      // (since bash throws error for unprotected parentheses)
      if((summons[0] == '\"' && summons[summons.size()-1] == '\"') || (summons[0] == '\'' && summons[summons.size()-1] == '\'')){
        summons.erase(summons.begin()); summons.erase(summons.begin()+summons.size()-1);
      }
    }
    return AFLUXCall(summons);
  }
}

namespace aflowlib {
  string AFLUXCall(vector<string>& matchbook){

    // Performs AFLUX call based on vector of matchbook entries
    string summons = aurostd::joinWDelimiter(matchbook,",");

    return AFLUXCall(summons);
  }
}

namespace aflowlib {
  string AFLUXCall(string& summons){

    // Performs AFLUX call based on summons input

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "AFLUXCall()";

    // percent encoding (otherwise it will not work)
    // NOT NEEDED - aurostd::StringSubst(summons,"\'","%27"); // percent encoding for "'" 
    aurostd::StringSubst(summons," ","%20");  // percent encoding for space 
    aurostd::StringSubst(summons,"#","%23");  // percent encoding for "#"

    if(LDEBUG) {
      cerr << function_name << ": Summons = " << summons << endl;
      cerr << function_name << ": Peforming call ... please be patient ..." << endl;
    }

    aflowlib::APIget API_socket(summons);
    stringstream response; response << API_socket;
    return response.str();

  }
}

namespace aflowlib {
  vector<vector<std::pair<string,string> > > getPropertiesFromAFLUXResponse(string& response){

    // Puts list of keyword-value pairs into a vector corresponding to each entry
    // Assumes the response format to be "format(aflow)", i.e., "|" delimiter
    // Here, pair.first=<keyword> and pair.second=<value>
    // In order to be general, all keywords and values are stored as a string 

    string function_name = XPID + "aflowlib::getPropertiesFromAFLUXResponse()";

    vector<vector<std::pair<string,string> > > properties_response;

    vector<string> entries,fields,key_value;
    aurostd::string2tokens(response,entries,"\n");

    // for each entry in response
    for(uint e=0;e<entries.size();e++){
      // split into key-value pairs
      aurostd::string2tokens(entries[e],fields,"|");

      // properties for a particular entry
      vector<std::pair<string,string> > property_pairs;
      for(uint i=0;i<fields.size();i++){
        aurostd::string2tokens(fields[i],key_value,"=");
        if(key_value.size()!=2 && !aurostd::substring2bool(fields[i],"example") && !aurostd::substring2bool(fields[i],"description")){ 
          cerr << function_name << "::ERROR: Cannot find key-value pair splitting on \"=\" for the following field: \"" << fields[i] << "\"." << endl; 
          exit(1);
        }
        std::pair<string,string> property; 
        property.first = aurostd::RemoveWhiteSpaces(key_value[0]);  // key 
        property.second = aurostd::RemoveWhiteSpaces(key_value[1]); // value
        property_pairs.push_back(property);
      }
      properties_response.push_back(property_pairs);
    }
    return properties_response;
  }
}

//DX+FR20190206 - AFLUX functionality via command line - END





// ***************************************************************************
//[SC20200327 - OBSOLETE]namespace aflowlib {
//[SC20200327 - OBSOLETE]  uint WEB_Aflowlib_Entry_PHP(string options,ostream& oss) {
//[SC20200327 - OBSOLETE]    bool LDEBUG=(FALSE || XHOST.DEBUG);
//[SC20200327 - OBSOLETE]    string soliloquy="aflowlib::WEB_Aflowlib_Entry():";
//[SC20200327 - OBSOLETE]    if(LDEBUG) cout << "aflowlib::WEB_Aflowlib_Entry: begin<br>" << endl;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]    stringstream num_prec;
//[SC20200327 - OBSOLETE]    vector<string> voptions;
//[SC20200327 - OBSOLETE]    aurostd::string2tokens(options,voptions,",");
//[SC20200327 - OBSOLETE]    if(voptions.size()==0) {
//[SC20200327 - OBSOLETE]      init::ErrorOption(options,"aflowlib::WEB_Aflowlib_Entry","aflow --aflowlib=entry");
//[SC20200327 - OBSOLETE]      exit(0);
//[SC20200327 - OBSOLETE]    } 
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]    // move on
//[SC20200327 - OBSOLETE]    for(uint ioption=0;ioption<voptions.size();ioption++) {
//[SC20200327 - OBSOLETE]      //  oss << voptions.at(ioption) << endl;
//[SC20200327 - OBSOLETE]      string option=voptions.at(ioption); //aurostd::args2attachedstring(argv,"--aflowlib=",(string) "nan");
//[SC20200327 - OBSOLETE]      if(option.at(option.size()-1)=='/'|| option.at(option.size()-1)=='.') option.erase(option.end()-1,option.end()-0); //  some demoronization
//[SC20200327 - OBSOLETE]      if(option.at(0)=='/'|| option.at(0)=='.') option.erase(option.begin(),option.begin()+1); //  some demoronization
//[SC20200327 - OBSOLETE]      string directory="";
//[SC20200327 - OBSOLETE]      string directory_RAW="",directory_LIB="",directory_WEB="";
//[SC20200327 - OBSOLETE]      string directory_AUID_LIB="",directory_AUID_RAW="",directory_AUID_WEB="";
//[SC20200327 - OBSOLETE]      string url_WEB;
//[SC20200327 - OBSOLETE]      string label="";
//[SC20200327 - OBSOLETE]      //string line_gif="<br><img border=0 width=60% height=2 src=http://materials.duke.edu/auro/images/line.gif><br><br>";
//[SC20200327 - OBSOLETE]      string line_rule="<hr width=\"60%\" style=\"background:black; border:0; height:2px; text-align:left;margin-left:0\" /><br>";
//[SC20200327 - OBSOLETE]      string art058_link=" [<a href=https://doi.org/10.1016/j.commatsci.2010.05.010 target=\"_blank\"><font color=black><i>cite</i></font></a>]";
//[SC20200327 - OBSOLETE]      string art064_link=" [<a href=https://doi.org/10.1021/co200012w target=\"_blank\"><font color=black><i>cite</i></font></a>]";
//[SC20200327 - OBSOLETE]      string icsd_link=" [<a href=https://www.fiz-karlsruhe.com/icsd.html target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string aflow_ael_readme=" [<a href=http://materials.duke.edu/AFLOW/README_AFLOW_AEL.TXT target=\"_blank\"><font color=black><i>info</i></font></a>]"; //CO20180817
//[SC20200327 - OBSOLETE]      string art096_link=" [<a href=https://doi.org/10.1103/PhysRevB.90.174107 target=\"_blank\"><font color=black><i>cite</i></font></a>]";
//[SC20200327 - OBSOLETE]      string art100_link=" [<a href=https://www.nature.com/articles/sdata20159 target=\"_blank\"><font color=black><i>cite</i></font></a>]";
//[SC20200327 - OBSOLETE]      string aflow_agl_readme=" [<a href=http://materials.duke.edu/AFLOW/README_AFLOW_AGL.TXT target=\"_blank\"><font color=black><i>info</i></font></a>]"; //CO20180817
//[SC20200327 - OBSOLETE]      string art115_link=" [<a href=https://doi.org/10.1103/PhysRevMaterials.1.015401 target=\"_blank\"><font color=black><i>cite</i></font></a>]"; //CO20180817
//[SC20200327 - OBSOLETE]      string aflow_sym_readme=" [<a href=http://materials.duke.edu/AFLOW/README_AFLOW_SYM.TXT target=\"_blank\"><font color=black><i>info</i></font></a>]"; //CO20180817
//[SC20200327 - OBSOLETE]      string art135_link=" [<a href=https://doi.org/10.1107/S2053273318003066 target=\"_blank\"><font color=black><i>cite</i></font></a>]"; //CO20180817
//[SC20200327 - OBSOLETE]      int atomCOUNT=0;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      //DX20180817
//[SC20200327 - OBSOLETE]      string bravais_lattice_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string bravais_lattice_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string lattice_system_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string lattice_variation_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string lattice_system_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string lattice_variation_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string Pearson_symbol_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string Pearson_symbol_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string sg_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string sg2_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg2 target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string spacegroup_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]      string spacegroup_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      aflowlib::_aflowlib_entry aentry;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      xoption vflags;
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::PREAMBLE",TRUE);
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::CALCULATION",TRUE);
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::JMOL",TRUE);
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::EDATA_ORIG",FALSE);
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::EDATA_RELAX",TRUE);
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::THERMODYNAMICS",TRUE);
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::MAGNETIC",TRUE);
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::ELECTRONIC",FALSE);     // will setup later
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::SCINTILLATION",TRUE);   // will setup later
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::AGL",FALSE);            // will setup later
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::AEL",FALSE);            // will setup later
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::BADER",FALSE);          // will setup later
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // check if ICSD inside (anyway)
//[SC20200327 - OBSOLETE]      string lattices[]={"BCC","BCT","CUB","FCC","HEX","MCL","MCLC","ORC","ORCC","ORCF","ORCI","RHL","TET","TRI"};
//[SC20200327 - OBSOLETE]      vector<string> vline,tokens;
//[SC20200327 - OBSOLETE]      // [OBSOLETE]   vector<string> tokens;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      string html_TAB=" target=\"_blank\"";
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: [1]<br>" << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: [4]<br>" << endl;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: option=" << option << endl;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // START SEARCH
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::FOUND",FALSE);
//[SC20200327 - OBSOLETE]      string catalog="",auid="";
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // *****************************************************
//[SC20200327 - OBSOLETE]      if(!vflags.flag("FLAG::FOUND") && aurostd::substring2bool(aurostd::tolower(option),"aflow:")) { // CHECK AUID
//[SC20200327 - OBSOLETE]        if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: option=" << option << endl;
//[SC20200327 - OBSOLETE]        string auid=aurostd::tolower(option);
//[SC20200327 - OBSOLETE]        if(auid.size()!=22) {
//[SC20200327 - OBSOLETE]          cerr << "WEB_Aflowlib_Entry_PHP: error on size of auid=" << auid << endl;
//[SC20200327 - OBSOLETE]          oss << "WEB_Aflowlib_Entry_PHP: error on size of auid=" << auid << endl;
//[SC20200327 - OBSOLETE]          exit(0);
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        // OLD
//[SC20200327 - OBSOLETE]        // directory_AUID_LIB=init::AFLOW_Projects_Directories("AUID")+"/LIB/"+auid.substr(0,8); for(uint i=8;i<=20;i+=2) directory_AUID_LIB+="/"+auid.substr(i,2);  // splitting aflow:ab/cd..
//[SC20200327 - OBSOLETE]        // directory_AUID_RAW=init::AFLOW_Projects_Directories("AUID")+"/RAW/"+auid.substr(0,8); for(uint i=8;i<=20;i+=2) directory_AUID_RAW+="/"+auid.substr(i,2);  // splitting aflow:ab/cd..
//[SC20200327 - OBSOLETE]        // directory_AUID_WEB=init::AFLOW_Projects_Directories("AUID")+"/WEB/"+auid.substr(0,8); for(uint i=8;i<=20;i+=2) directory_AUID_WEB+="/"+auid.substr(i,2);  // splitting aflow:ab/cd..
//[SC20200327 - OBSOLETE]        // NEW
//[SC20200327 - OBSOLETE]        directory_AUID_LIB=init::AFLOW_Projects_Directories("AUID")+"/"+auid.substr(0,8); for(uint i=8;i<=20;i+=2) directory_AUID_LIB+="/"+auid.substr(i,2);  // splitting aflow:ab/cd..
//[SC20200327 - OBSOLETE]        directory_AUID_WEB=directory_AUID_LIB+"/WEB";
//[SC20200327 - OBSOLETE]        directory_AUID_RAW=directory_AUID_LIB+"/RAW";
//[SC20200327 - OBSOLETE]        directory_AUID_LIB=directory_AUID_LIB+"/LIB";
//[SC20200327 - OBSOLETE]        if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: directory_AUID_LIB=" << directory_AUID_LIB << endl;
//[SC20200327 - OBSOLETE]        if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: directory_AUID_RAW=" << directory_AUID_RAW << endl;
//[SC20200327 - OBSOLETE]        if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: directory_AUID_WEB=" << directory_AUID_WEB << endl;
//[SC20200327 - OBSOLETE]        directory="";
//[SC20200327 - OBSOLETE]        if(!aurostd::FileExist(directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]          cerr << "WEB_Aflowlib_Entry_PHP: entry does not exist =" << directory_AUID_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << endl;exit(0);
//[SC20200327 - OBSOLETE]          oss << "WEB_Aflowlib_Entry_PHP: entry does not exist =" << directory_AUID_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << endl;exit(0);
//[SC20200327 - OBSOLETE]        } else {
//[SC20200327 - OBSOLETE]          _aflowlib_entry entry_tmp(string(directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT));
//[SC20200327 - OBSOLETE]          directory=entry_tmp.aurl;
//[SC20200327 - OBSOLETE]          auid=entry_tmp.aurl;
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"aflowlib.duke.edu:","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"materials.duke.edu:","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/ICSD_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/ICSD_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB0_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB0_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB1_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB1_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB2_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB2_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB3_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB3_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB4_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB4_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB5_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB5_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB6_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB6_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB7_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB7_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB8_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB8_WEB/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/LIB9_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB9_WEB/","");
//[SC20200327 - OBSOLETE]          vflags.flag("FLAG::AUID",TRUE);
//[SC20200327 - OBSOLETE]          vflags.flag("FLAG::FOUND",TRUE);
//[SC20200327 - OBSOLETE]          catalog=entry_tmp.catalog;
//[SC20200327 - OBSOLETE]          label=directory;
//[SC20200327 - OBSOLETE]          for(uint ilat=0;ilat<14;ilat++)
//[SC20200327 - OBSOLETE]            aurostd::StringSubst(label,lattices[ilat]+"/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(label,"/",".");
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // *****************************************************
//[SC20200327 - OBSOLETE]      if(!vflags.flag("FLAG::FOUND")) { // tests with proto name
//[SC20200327 - OBSOLETE]        string dir2test;
//[SC20200327 - OBSOLETE]        dir2test=option;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        vector<string> vdir2test;
//[SC20200327 - OBSOLETE]        for(uint i0=0;i0<dir2test.size();i0++) { // allow all identical indices...
//[SC20200327 - OBSOLETE]          for(uint i1=i0;i1<dir2test.size();i1++) { // allow all identical indices
//[SC20200327 - OBSOLETE]            for(uint i2=i1;i2<dir2test.size();i2++) { // allow all identical indices
//[SC20200327 - OBSOLETE]              //	    for(uint i3=i2;i3<dir2test.size();i3++) { // allow all identical indices
//[SC20200327 - OBSOLETE]              dir2test=option;
//[SC20200327 - OBSOLETE]              if(dir2test.at(i0)=='.') dir2test.at(i0)='/';
//[SC20200327 - OBSOLETE]              if(dir2test.at(i1)=='.') dir2test.at(i1)='/';
//[SC20200327 - OBSOLETE]              if(dir2test.at(i2)=='.') dir2test.at(i2)='/';
//[SC20200327 - OBSOLETE]              //      if(dir2test.at(i3)=='.') dir2test.at(i3)='/';
//[SC20200327 - OBSOLETE]              bool found=false;
//[SC20200327 - OBSOLETE]              for(uint i=0;i<vdir2test.size()&&!found;i++) if(dir2test==vdir2test.at(i)) found=true;
//[SC20200327 - OBSOLETE]              if(!found) {
//[SC20200327 - OBSOLETE]                vdir2test.push_back(dir2test);
//[SC20200327 - OBSOLETE]                //		cerr << vdir2test.size() << endl;
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]            //	  }
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: testing dir2test=" << dir2test << endl;
//[SC20200327 - OBSOLETE]        for(uint i=0;i<vdir2test.size()&&!vflags.flag("FLAG::FOUND");i++) { // allow i=j=k so that 1 OR 2 OR 3 dots are also tested
//[SC20200327 - OBSOLETE]          dir2test=vdir2test.at(i);
//[SC20200327 - OBSOLETE]          //		cout << "testing(" << i << ") = " << dir2test << endl;
//[SC20200327 - OBSOLETE]          for(uint ilat=0;ilat<14&&!vflags.flag("FLAG::FOUND");ilat++) {
//[SC20200327 - OBSOLETE]            //	cerr << string(init::AFLOW_Projects_Directories("ICSD")+"/RAW/"+lattices[ilat]+"/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT) << endl;
//[SC20200327 - OBSOLETE]            if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("ICSD")+"/RAW/"+lattices[ilat]+"/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]              catalog="ICSD";directory=lattices[ilat]+"/"+dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB0")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB0";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB1")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB1";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB2";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB3")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB3";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB4")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB4";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB5")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB5";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB6")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB6";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB7")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB7";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB8")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB8";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB9")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]            catalog="LIB9";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // **********************************************************************************************************
//[SC20200327 - OBSOLETE]      // TRY ICSD LINK
//[SC20200327 - OBSOLETE]      // **********************************************************************************************************
//[SC20200327 - OBSOLETE]      if(!vflags.flag("FLAG::FOUND")) { // icsd link
//[SC20200327 - OBSOLETE]        string directory_ICSD2LINK=init::AFLOW_Projects_Directories("AUID")+"/icsd:/"+option;
//[SC20200327 - OBSOLETE]        aurostd::StringSubst(directory_ICSD2LINK,"ICSD:","icsd:");
//[SC20200327 - OBSOLETE]        aurostd::StringSubst(directory_ICSD2LINK,"icsd:icsd:","icsd:");    
//[SC20200327 - OBSOLETE]        //	cerr << directory_ICSD2LINK << endl;
//[SC20200327 - OBSOLETE]        if(aurostd::FileExist(directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
//[SC20200327 - OBSOLETE]          _aflowlib_entry entry_tmp(string(directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT));
//[SC20200327 - OBSOLETE]          directory=entry_tmp.aurl;
//[SC20200327 - OBSOLETE]          auid=entry_tmp.aurl;
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"aflowlib.duke.edu:","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"materials.duke.edu:","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(directory,"AFLOWDATA/ICSD_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/ICSD_WEB/","");
//[SC20200327 - OBSOLETE]          vflags.flag("FLAG::ICSD",TRUE);
//[SC20200327 - OBSOLETE]          vflags.flag("FLAG::FOUND",TRUE);
//[SC20200327 - OBSOLETE]          catalog=entry_tmp.catalog;
//[SC20200327 - OBSOLETE]          label=directory;
//[SC20200327 - OBSOLETE]          for(uint ilat=0;ilat<14;ilat++)
//[SC20200327 - OBSOLETE]            aurostd::StringSubst(label,lattices[ilat]+"/","");
//[SC20200327 - OBSOLETE]          aurostd::StringSubst(label,"/",".");
//[SC20200327 - OBSOLETE]          //	  cerr << directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << endl;
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // **********************************************************************************************************
//[SC20200327 - OBSOLETE]      // SHOULD BE FOUND
//[SC20200327 - OBSOLETE]      // **********************************************************************************************************
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(catalog=="ICSD") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::ICSD",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("ICSD")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("ICSD")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("ICSD")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/ICSD_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(catalog=="LIB0") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB0",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB0")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB0")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB0")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB0_RAW/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB1") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB1",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB1")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB1")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB1")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB1_RAW/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB2") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB2",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB2")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+directory; // June 2016
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB2_RAW/"+directory; // May 2014
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB3") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB3",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB3")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB3")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB3")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB3_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB4") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB4",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB4")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB4")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB4")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB4_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB5") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB5",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB5")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB5")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB5")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB5_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB6") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB6",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB6")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB6")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB6")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB6_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB7") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB7",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB7")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB7")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB7")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB7_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB8") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB8",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB8")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB8")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB8")+"/RAW/"+directory;
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB8_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(catalog=="LIB9") {
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::LIB9",TRUE);
//[SC20200327 - OBSOLETE]        directory_LIB=init::AFLOW_Projects_Directories("LIB9")+"/LIB/"+directory;
//[SC20200327 - OBSOLETE]        directory_WEB=init::AFLOW_Projects_Directories("LIB9")+"/WEB/"+directory;
//[SC20200327 - OBSOLETE]        directory_RAW=init::AFLOW_Projects_Directories("LIB9")+"/RAW/"+directory;	
//[SC20200327 - OBSOLETE]        url_WEB="/AFLOWDATA/LIB9_WEB/"+directory;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // LDEBUG=1;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::FOUND\")=" << vflags.flag("FLAG::FOUND") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::ICSD\")=" << vflags.flag("FLAG::ICSD") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB1\")=" << vflags.flag("FLAG::LIB1") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB2\")=" << vflags.flag("FLAG::LIB2") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB3\")=" << vflags.flag("FLAG::LIB3") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB4\")=" << vflags.flag("FLAG::LIB4") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB5\")=" << vflags.flag("FLAG::LIB5") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB7\")=" << vflags.flag("FLAG::LIB6") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB7\")=" << vflags.flag("FLAG::LIB7") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB8\")=" << vflags.flag("FLAG::LIB8") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::LIB9\")=" << vflags.flag("FLAG::LIB9") << endl;
//[SC20200327 - OBSOLETE]      if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vflags.flag(\"FLAG::AUID\")=" << vflags.flag("FLAG::AUID") << endl;   
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // now start
//[SC20200327 - OBSOLETE]      // got it  ?
//[SC20200327 - OBSOLETE]      if(!aurostd::FileExist(directory_RAW+"/"+_AFLOWIN_)) directory_RAW="";
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(!directory.empty()) { // play with aentry.entry
//[SC20200327 - OBSOLETE]        aurostd::StringSubst(label,"/",".");
//[SC20200327 - OBSOLETE]        aentry.file2aflowlib(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,oss);  //   oss << aentry.entry << endl;
//[SC20200327 - OBSOLETE]        directory_AUID_LIB=init::AFLOW_Projects_Directories("AUID")+"/"+aflowlib::auid2directory(aentry.auid);
//[SC20200327 - OBSOLETE]        directory_AUID_WEB=directory_AUID_LIB+"/WEB";
//[SC20200327 - OBSOLETE]        directory_AUID_RAW=directory_AUID_LIB+"/RAW";
//[SC20200327 - OBSOLETE]        directory_AUID_LIB=directory_AUID_LIB+"/LIB";    
//[SC20200327 - OBSOLETE]        aurostd::string2tokens(aentry.sg2,tokens,"#");if(tokens.size()>0) aentry.sg2=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]        if(aentry.vfiles_WEB.size()==0) aentry.vfiles_WEB=aentry.vfiles;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // check AGL/AEL
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::ELECTRONIC",aurostd::substring2bool(aentry.vloop,"bands"));
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::SCINTILLATION",aurostd::substring2bool(aentry.vloop,"bands"));
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::AGL",aurostd::substring2bool(aentry.vloop,"agl"));
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::AEL",aurostd::substring2bool(aentry.vloop,"ael"));
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::BADER",aurostd::substring2bool(aentry.vloop,"bader"));
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(XHOST.hostname=="nietzsche.mems.duke.edu") {
//[SC20200327 - OBSOLETE]        oss << "<b>DEBUG: only in " << XHOST.hostname << "</b><br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "XHOST.hostname=" << XHOST.hostname << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "auid=" << aentry.auid << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "option=" << option << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "label=" << label << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "directory=" << directory << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "catalog=" << aentry.catalog << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::ICSD\")=" << vflags.flag("FLAG::ICSD") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB0\")=" << vflags.flag("FLAG::LIB0") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB1\")=" << vflags.flag("FLAG::LIB1") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB2\")=" << vflags.flag("FLAG::LIB2") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB3\")=" << vflags.flag("FLAG::LIB3") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB4\")=" << vflags.flag("FLAG::LIB4") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB5\")=" << vflags.flag("FLAG::LIB5") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB6\")=" << vflags.flag("FLAG::LIB6") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB7\")=" << vflags.flag("FLAG::LIB7") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB8\")=" << vflags.flag("FLAG::LIB8") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::LIB9\")=" << vflags.flag("FLAG::LIB9") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::AUID\")=" << vflags.flag("FLAG::AUID") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::FOUND\")=" << vflags.flag("FLAG::FOUND") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "directory_LIB=" << directory_LIB << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "directory_WEB=" << directory_WEB << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "directory_RAW=" << directory_RAW << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "directory_AUID_LIB=" << directory_AUID_LIB << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "directory_AUID_WEB=" << directory_AUID_WEB << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "directory_AUID_RAW=" << directory_AUID_RAW << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::PREAMBLE\")=" << vflags.flag("FLAG::PREAMBLE")  << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::CALCULATION\")=" << vflags.flag("FLAG::CALCULATION") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::JMOL\")=" << vflags.flag("FLAG::JMOL")  << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::EDATA_ORIG\")=" << vflags.flag("FLAG::EDATA_ORIG") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::EDATA_RELAX\")=" << vflags.flag("FLAG::EDATA_RELAX") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::THERMODYNAMICS\")=" << vflags.flag("FLAG::THERMODYNAMICS") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::MAGNETIC\")=" << vflags.flag("FLAG::MAGNETIC") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::ELECTRONIC\")=" << vflags.flag("FLAG::ELECTRONIC") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::SCINTILLATION\")=" << vflags.flag("FLAG::SCINTILLATION") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::AGL\")=" << vflags.flag("FLAG::AGL") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::AEL\")=" << vflags.flag("FLAG::AEL") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "vflags.flag(\"FLAG::BADER\")=" << vflags.flag("FLAG::BADER") << "<br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "aentry.loop=" << aentry.loop << "<br>" << endl;
//[SC20200327 - OBSOLETE]        // oss << "aflowlib.out=" << aurostd::efile2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT) << "<br>" << endl;
//[SC20200327 - OBSOLETE]        // oss << "aflowlib.json=" << aurostd::efile2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON) << "<br>" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      //CO20180523 - fixing for LIB6 missing from /www directory
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(XHOST.hostname, "aflowlib")){
//[SC20200327 - OBSOLETE]        vflags.flag("FLAG::JMOL",FALSE);
//[SC20200327 - OBSOLETE]        string web_path="/www"+url_WEB;
//[SC20200327 - OBSOLETE]        if(LDEBUG) {cerr << soliloquy << " web_path=" << web_path << endl;}
//[SC20200327 - OBSOLETE]        if(aurostd::IsDirectory(web_path)){
//[SC20200327 - OBSOLETE]          vector<string> web_path_vfiles;
//[SC20200327 - OBSOLETE]          aurostd::DirectoryLS(web_path,web_path_vfiles);
//[SC20200327 - OBSOLETE]          if(aurostd::substring2bool(web_path_vfiles,label+".cif")){vflags.flag("FLAG::JMOL",TRUE);}
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        if(LDEBUG) {cerr << soliloquy << " vflags.flag(\"FLAG::JMOL\")=" << vflags.flag("FLAG::JMOL") << endl;}
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // make ORIG vs RELAX
//[SC20200327 - OBSOLETE]      vflags.flag("FLAG::EDATA_ORIG",FALSE);
//[SC20200327 - OBSOLETE]      if(aentry.Bravais_lattice_orig!=aentry.Bravais_lattice_relax && aentry.Bravais_lattice_orig!="nan" && aentry.Bravais_lattice_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);
//[SC20200327 - OBSOLETE]      if(aentry.lattice_variation_orig!=aentry.lattice_variation_relax && aentry.lattice_variation_orig!="nan" && aentry.lattice_variation_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);
//[SC20200327 - OBSOLETE]      if(aentry.lattice_system_orig!=aentry.lattice_system_relax && aentry.lattice_system_orig!="nan" && aentry.lattice_system_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);
//[SC20200327 - OBSOLETE]      if(aentry.Pearson_symbol_orig!=aentry.Pearson_symbol_relax && aentry.Pearson_symbol_orig!="nan" && aentry.Pearson_symbol_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // not found
//[SC20200327 - OBSOLETE]      // oss << "Thank you for your query, but unfortunately entry " << label << " has not yet been calculated. If you want to report this omission please email Dr. Rose at aflowdev@aflowlib.duke.ed.<br>" << endl;
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // PREAMBLE BEGIN
//[SC20200327 - OBSOLETE]      string title=label;
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::PREAMBLE") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        if(vflags.flag("FLAG::ICSD")) {
//[SC20200327 - OBSOLETE]          aurostd::string2tokens(label,tokens,"_");
//[SC20200327 - OBSOLETE]          title=tokens.at(0);
//[SC20200327 - OBSOLETE]          for(uint i=0;i<=9;i++) aurostd::StringSubst(title,aurostd::utype2string<uint>(i),string("<sub>"+aurostd::utype2string<uint>(i)+"</sub>"));
//[SC20200327 - OBSOLETE]          title+=" (ICSD# "+tokens.at(2)+")";//+directory;
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::ELECTRONIC")){ //CO20180502
//[SC20200327 - OBSOLETE]        oss.setf(std::ios::fixed,std::ios::floatfield);
//[SC20200327 - OBSOLETE]        oss.precision(4);
//[SC20200327 - OBSOLETE]        oss << "<script type=\"text/javascript\">" << endl;
//[SC20200327 - OBSOLETE]        aurostd::xoption aaa;
//[SC20200327 - OBSOLETE]        stringstream bandsdata;
//[SC20200327 - OBSOLETE]        oss << "var d3_bands_data = "; estructure::BANDSDATA_JSON(aaa, directory_LIB, bandsdata,true); oss << bandsdata.str();  //GG
//[SC20200327 - OBSOLETE]        oss << ";</script>" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      oss << "<! HS WORK BEFORE HERE> " << endl;
//[SC20200327 - OBSOLETE]      oss << "<div id=\"content\">" << endl;
//[SC20200327 - OBSOLETE]      oss << "<div class=\"title\">" << endl;
//[SC20200327 - OBSOLETE]      oss << "<FONT SIZE=+3> " << title << " </FONT></div>" << endl;
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // COMPOUND
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- compound: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<FONT SIZE=+3> " << aentry.compound << " </FONT><br>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- compound: END -->" << endl;
//[SC20200327 - OBSOLETE]      }	
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // LICENSE
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        string LICENSE="The data included within the aflow.org repository is free for scientific, academic and non-commercial purposes. Any other use is prohibited.";
//[SC20200327 - OBSOLETE]        oss << "<!-- license: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<div class = \"url_text\">" << endl;
//[SC20200327 - OBSOLETE]        oss << "<br><b>LICENSE: <span class=\"url_text\">" << LICENSE << "</span></b>" << endl;
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- license: END -->" << endl;
//[SC20200327 - OBSOLETE]      }	
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // PREAMBLE BEGIN
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::PREAMBLE") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        // string URL=string("http://aflow.org/material.php?proto_name=")+label;
//[SC20200327 - OBSOLETE]        // string URL=string("http://aflow.org/material.php?id=")+label;   // with label
//[SC20200327 - OBSOLETE]        string URL=string("http://aflow.org/material.php?id=")+aentry.auid; // with auid
//[SC20200327 - OBSOLETE]        oss << "<!-- preamble: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        //      oss << "<FONT SIZE=+3> " << aentry.compound << " </FONT>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<div class = \"url\">" << endl;
//[SC20200327 - OBSOLETE]        oss << "Permanent URL: <span class=\"url_text\">" << URL << "</span>" << endl;
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<br><FONT SIZE=+0><b>aflow.org</b> web entry generator V" << string(AFLOW_VERSION) << " [built="  << TODAY << "]" << "</font>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- preamble: END -->" << endl;
//[SC20200327 - OBSOLETE]      }	
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // NEW JSMOL BEGIN
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::JMOL") && !directory.empty() && aurostd::substring2bool(aentry.vfiles_WEB,label+".cif")) {  //CO20180523 - fixing for LIB6 missing from /www directory
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]if label+".cif" is available, assume "_sprim" and "_sconv" are too
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]space group stuff found in bob's file now
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]oss << "<br><br><b>Space Group</b>: " << (aurostd::string2utype<int>(aentry.spacegroup_relax)>0 ? GetSpaceGroupName(aurostd::string2utype<int>(aentry.spacegroup_relax)) : "N/A" ) << "  (#" << aentry.spacegroup_relax << ")" << endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]if(aurostd::string2utype<int>(aentry.spacegroup_relax)>0) {
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]  oss << "<br><br><b>Space Group</b>: " <<  GetSpaceGroupName(aurostd::string2utype<int>(aentry.spacegroup_relax)) << "  (#" << aentry.spacegroup_relax << ")" << endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]} else {
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]  oss << "<br><br><b>Space Group</b>: " <<  "N/A" << "  (#" << aentry.spacegroup_relax << ")" << endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]}
//[SC20200327 - OBSOLETE]        oss << "<!-- jmol: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!--div class = \"jmol\"-->" << endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]oss << "<script type=\"text/javascript\" src=\"/search/Lib/JS/JSmol.min.js\"></script>" << endl;  //CO20170622
//[SC20200327 - OBSOLETE]        //[OBSOLETE CO20170628 - per Bob/JMOL]string JMOL_PATH="http://aflowlib.duke.edu/users/jmolers/test/jsmol";
//[SC20200327 - OBSOLETE]        string JMOL_PATH="http://aflowlib.duke.edu/search/Lib/JS/JSMol";
//[SC20200327 - OBSOLETE]        oss << "<script type=\"text/javascript\" src=\"" << JMOL_PATH << "/JSmol.min.js\"></script>" << endl;  //CO20170622
//[SC20200327 - OBSOLETE]        oss << "<script type=\"text/javascript\">" << endl;
//[SC20200327 - OBSOLETE]        //CO20170622 START
//[SC20200327 - OBSOLETE]        //build our standard AFLOW object, add from aentry as needed
//[SC20200327 - OBSOLETE]        oss << "AFLOW={};" << endl;
//[SC20200327 - OBSOLETE]        oss << "AFLOW.version = \"" << string(AFLOW_VERSION) << "\";" << endl;
//[SC20200327 - OBSOLETE]        oss << "AFLOW.url_WEB = \"http://aflowlib.duke.edu" << url_WEB << "\";" << endl; //CO check with bob, maybe aflowlib.duke.edu?
//[SC20200327 - OBSOLETE]        string system_name=KBIN::ExtractSystemName(directory_LIB);
//[SC20200327 - OBSOLETE]        //    cerr << system_name << endl; exit(0);
//[SC20200327 - OBSOLETE]        //  system_name="xxxx";
//[SC20200327 - OBSOLETE]        oss << "AFLOW.label = \"" << system_name << "\";" << endl;
//[SC20200327 - OBSOLETE]        //  oss << "AFLOW.label = \"" << "xxx" << "\";" << endl;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        oss << "AFLOW.spaceGroupNo = " << aentry.spacegroup_relax << ";" << endl;
//[SC20200327 - OBSOLETE]        oss << "AFLOW.spaceGroupName = \"" << GetSpaceGroupName(aurostd::string2utype<int>(aentry.spacegroup_relax)) << "\";" << endl;
//[SC20200327 - OBSOLETE]        double aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF;
//[SC20200327 - OBSOLETE]        cif2data(directory_WEB+"/"+label+"_sconv.cif",aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
//[SC20200327 - OBSOLETE]        oss << "AFLOW.cif_sconv = [" << aCIF << "," << bCIF << "," << cCIF << "," << alphaCIF << "," << betaCIF << "," << gammaCIF << "];" << endl;
//[SC20200327 - OBSOLETE]        cif2data(directory_WEB+"/"+label+".cif",aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
//[SC20200327 - OBSOLETE]        oss << "AFLOW.cif = [" << aCIF << "," << bCIF << "," << cCIF << "," << alphaCIF << "," << betaCIF << "," << gammaCIF << "];" << endl;
//[SC20200327 - OBSOLETE]        cif2data(directory_WEB+"/"+label+"_sprim.cif",aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
//[SC20200327 - OBSOLETE]        oss << "AFLOW.cif_sprim = [" << aCIF << "," << bCIF << "," << cCIF << "," << alphaCIF << "," << betaCIF << "," << gammaCIF << "];" << endl;
//[SC20200327 - OBSOLETE]        string sym2json; //PC+DX20180723
//[SC20200327 - OBSOLETE]        if(aurostd::FileExist(directory_RAW+"/"+"aflow.fgroup.bands.json.xz")) //PC+DX20180723
//[SC20200327 - OBSOLETE]        { aurostd::xzfile2string(directory_RAW+"/"+"aflow.fgroup.bands.json.xz",sym2json); //PC+DX20180723
//[SC20200327 - OBSOLETE]          oss << "AFLOW.sym2json ="; //PC+DX20180723
//[SC20200327 - OBSOLETE]          oss << sym2json;    //PC+DX20180723
//[SC20200327 - OBSOLETE]          oss << ";" << endl; } //PC+DX20180723
//[SC20200327 - OBSOLETE]        else if (aurostd::FileExist(directory_RAW+"/"+"aflow.fgroup.relax.json.xz")) //PC+DX20180723
//[SC20200327 - OBSOLETE]        { aurostd::xzfile2string(directory_RAW+"/"+"aflow.fgroup.relax.json.xz",sym2json); //PC+DX20180723
//[SC20200327 - OBSOLETE]          oss << "AFLOW.sym2json ="; //PC+DX20180723
//[SC20200327 - OBSOLETE]          oss << sym2json;    //PC+DX20180723
//[SC20200327 - OBSOLETE]          oss << ";" << endl; }//PC+DX20180723
//[SC20200327 - OBSOLETE]        else {cerr << "error" << endl; //PC+DX20180723
//[SC20200327 - OBSOLETE]        }; //PC+DX20180723
//[SC20200327 - OBSOLETE]        //BEGIN BADER ISOSURFACES
//[SC20200327 - OBSOLETE]        if(vflags.flag("FLAG::BADER")){ //did we calculate bader?
//[SC20200327 - OBSOLETE]          if(aurostd::substring2bool(aentry.vfiles_WEB,label+"_Bader_20_"+aentry.vspecies.at(0)+".jvxl")) { //quick (not robust) test that bader loop ran fine
//[SC20200327 - OBSOLETE]            if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax")) {  //check that we have the right structure
//[SC20200327 - OBSOLETE]              xstructure xstr(directory_WEB+"/CONTCAR.relax",IOAFLOW_AUTO);
//[SC20200327 - OBSOLETE]              xstr.ReScale(1.0);
//[SC20200327 - OBSOLETE]              oss << "AFLOW.baderUnitcell = [";
//[SC20200327 - OBSOLETE]              string sep = "";
//[SC20200327 - OBSOLETE]              for(uint i=1;i<=3;i++)
//[SC20200327 - OBSOLETE]                for(uint j=1;j<=3;j++) {
//[SC20200327 - OBSOLETE]                  oss << sep << xstr.lattice(i,j);
//[SC20200327 - OBSOLETE]                  sep = ",";
//[SC20200327 - OBSOLETE]                }
//[SC20200327 - OBSOLETE]              oss << "];" << endl;
//[SC20200327 - OBSOLETE]              oss << "AFLOW.baderVSpecies = [" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(aentry.vspecies,"\""),",") << "];" << endl;
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        //END BADER ISOSURFACES
//[SC20200327 - OBSOLETE]        oss << "AFLOW.jsmolDir = \"" << JMOL_PATH << "\";" << endl;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        //adding bob's stuff
//[SC20200327 - OBSOLETE]        //string aflow_entry_js=AFLOW_ENTRY_JS;
//[SC20200327 - OBSOLETE]        oss << AFLOW_WEBAPP_ENTRY_JS;  //CO20170622 
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        oss << "</script>" << endl;
//[SC20200327 - OBSOLETE]        oss << "</td></tr></table>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!--/div-->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- jmol: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ELECTRONIC") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- GG bands: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<script type=\"text/javascript\" src=\"./Lib/JS/d3.min.js\"></script>" << endl;  //CO20170622  ///www/search/Lib/JS/d3.min.js
//[SC20200327 - OBSOLETE]        oss << "<script type=\"text/javascript\">" << endl;
//[SC20200327 - OBSOLETE]        oss << AFLOW_WEBAPP_BANDS_JS; //PC20180515
//[SC20200327 - OBSOLETE]        oss << "</script>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- GG bands: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // NEW JSMOL END
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // CALCULATION
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::CALCULATION") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Calculation properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\">Calculation details</h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //      oss << "<div class=\"calculation_details\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.code.empty()) {
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ":</h5><span class=\"value\">" << "[" << "<a href=\"" << url_WEB << "/\"" << html_TAB << ">entry</a>|" << "<a href=\"" << url_WEB << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "\"" << html_TAB << ">raw</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << ":</h5><span class=\"value\">" << "[" << "<a href=\"" << url_WEB << "/?format=json\"" << html_TAB << ">entry</a>|" << "<a href=\"" << url_WEB << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "\"" << html_TAB << ">raw</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        if(!aentry.auid.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW-UID [<a href=\"" << url_WEB << "/?auid\">auid</a>]:</h5><span class=\"value\">" << aentry.auid << "</span></div></div>" << endl; //PC20180515  //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.aurl.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW-URL [<a href=\"" << url_WEB << "/?aurl\">aurl</a>]:</h5><span class=\"value\">" << aentry.aurl << "</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.code.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"><i>Ab-initio</i> code [<a href=\"" << url_WEB << "/?code\">code</a>]:</h5><span class=\"value\">" << aentry.code << (aentry.calculation_cores>1?" (MPI) ":"") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.compound.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Composition per cell [<a href=\"" << url_WEB << "/?compound\">compound</a>]:</h5><span class=\"value\">" << aentry.compound << "</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.icsd.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">ICSD entry" << icsd_link << ":</h5><span class=\"value\">" << aentry.icsd << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.dft_type.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Pseudopotentials type [<a href=\"" << url_WEB << "/?pp_type\">pp_type</a>]:</h5><span class=\"value\">" << aentry.dft_type << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        for(uint i=0;i<aentry.vspecies_pp.size();i++)  
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">PP - Species&middot;Version [<a href=\"" << url_WEB << "/?species_pp_version\">species_pp_version</a>]:</h5><span class=\"value\">" << aentry.vspecies_pp_version.at(i) << "</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        for(uint i=0;i<aentry.vspecies_pp_ZVAL.size();i++)  
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">PP - Species&middot;ZVAL [<a href=\"" << url_WEB << "/?species_pp_ZVAL\">species_pp_ZVAL</a>]:</h5><span class=\"value\">" << aentry.vspecies_pp_ZVAL.at(i) << "</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">LDAU [T;{L};{U};{J}] [<a href=\"" << url_WEB << "/?ldau_TLUJ\">ldau_TLUJ</a>]:</h5><span class=\"value\">" << (!aentry.ldau_TLUJ.empty()?aentry.ldau_TLUJ+"  (type;{angular-1};{eV};{eV}) ":"no-LDAU") << "</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.calculation_time*aentry.calculation_cores>0.0) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Total CPU&middot;hours [<a href=\"" << url_WEB << "/?calculation_time\">calculation_time</a>*<a href=\"" << url_WEB << "/?calculation_cores\">calculation_cores</a>]/3600:</h5><span class=\"value\">" << aentry.calculation_time*aentry.calculation_cores/3600 << " hours</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.calculation_time>0.0) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Total Wall-time [<a href=\"" << url_WEB << "/?calculation_time\">calculation_time/3600</a>]:</h5><span class=\"value\">" << aentry.calculation_time/3600 << " hours</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.calculation_memory>0.0)
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Memory Used [<a href=\"" << url_WEB << "/?calculation_memory\">calculation_memory</a>]:</h5><span class=\"value\">" << aentry.calculation_memory << " MB</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.calculation_cores>0) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Number of cores [<a href=\"" << url_WEB << "/?calculation_cores\">calculation_cores</a>]:</h5><span class=\"value\">" << aentry.calculation_cores << (aentry.calculation_cores>1?" (MPI) ":"(serial)") << "</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.node_CPU_Model.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Processor - CPU Model [<a href=\"" << url_WEB << "/?node_CPU_Model\">node_CPU_Model</a>]:</h5><span class=\"value\">" << (!aentry.node_CPU_Model.empty()?aentry.node_CPU_Model:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.node_CPU_MHz>0.0) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Processor - CPU MHz [<a href=\"" << url_WEB << "/?node_CPU_MHz\">node_CPU_MHz</a>]:</h5><span class=\"value\">" << (aentry.node_CPU_MHz>1.0?aurostd::utype2string(aentry.node_CPU_MHz)+"MHz ":"unavailable") << " </span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.aflow_version.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?aflow_version\">aflow_version</a>]:</h5><span class=\"value\">" << (!aentry.aflow_version.empty()?aentry.aflow_version:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.catalog.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Catalog [<a href=\"" << url_WEB << "/?catalog\">catalog</a>]:</h5><span class=\"value\">" << (!aentry.catalog.empty()?aentry.catalog:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.data_api.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Data API [<a href=\"" << url_WEB << "/?data_api\">data_api</a>]:</h5><span class=\"value\">" << (!aentry.data_api.empty()?aentry.data_api:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.data_source.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Data Source [<a href=\"" << url_WEB << "/?data_source\">data_source</a>]:</h5><span class=\"value\">" << (!aentry.data_source.empty()?aentry.data_source:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.data_language.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Data Language [<a href=\"" << url_WEB << "/?data_language\">data_language</a>]:</h5><span class=\"value\">" << (!aentry.data_language.empty()?aentry.data_language:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.error_status.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Error Status [<a href=\"" << url_WEB << "/?error_status\">error_status</a>]:</h5><span class=\"value\">" << (!aentry.error_status.empty()?aentry.error_status:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.loop.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW loops [<a href=\"" << url_WEB << "/?loop\">loop</a>]:</h5><span class=\"value\">" << aentry.loop << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //     oss << "<li><span class=\"description\">AFLOW precision:</span>" << "unavailable" << "</li>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //    oss << "<li><span class=\"description\">AFLOW KPPRA:</span>" << "unavailable" << "</li>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.aflowlib_version.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOWLIB_entry version [<a href=\"" << url_WEB << "/?aflowlib_version\">aflowlib_version</a>]:</h5><span class=\"value\">" << (!aentry.aflowlib_version.empty()?aentry.aflowlib_version:"unavailable") << "</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.vaflowlib_date.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOWLIB_entry date [<a href=\"" << url_WEB << "/?aflowlib_date\">aflowlib_date</a>]:</h5><span class=\"value\">" << (!aentry.vaflowlib_date.empty()?aurostd::joinWDelimiter(aentry.vaflowlib_date,","):"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731  //CO20200624 - adding LOCK date
//[SC20200327 - OBSOLETE]        if(!aentry.author.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?author\">author</a>]:</h5><span class=\"value\">" << (!aentry.author.empty()?aentry.author:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.corresponding.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?corresponding\">corresponding</a>]:</h5><span class=\"value\">" << (!aentry.corresponding.empty()?aentry.corresponding:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(!aentry.sponsor.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?sponsor\">sponsor</a>]:</h5><span class=\"value\">" << (!aentry.sponsor.empty()?aentry.sponsor:"unavailable") << "</span></div></div>" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        // CORMAC
//[SC20200327 - OBSOLETE]        if(aentry.energy_cutoff!=AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">energy_cutoff [<a href=\"" << url_WEB << "/?energy_cutoff\">energy_cutoff</a>]:</h5><span class=\"value\">" << aentry.energy_cutoff << " eV</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        //      oss.precision(6);
//[SC20200327 - OBSOLETE]        if(aentry.delta_electronic_energy_convergence!=AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">delta_electronic_energy_convergence [<a href=\"" << url_WEB << "/?delta_electronic_energy_convergence\">delta_electronic_energy_convergence</a>]:</h5><span class=\"value\">" << 1000*aentry.delta_electronic_energy_convergence << " meV</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.delta_electronic_energy_threshold!=AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">delta_electronic_energy_threshold_ [<a href=\"" << url_WEB << "/?delta_electronic_energy_threshold\">delta_electronic_energy_threshold</a>]:</h5><span class=\"value\">" << 1000*aentry.delta_electronic_energy_threshold << " meV</span></div></div>" << endl; //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.nkpoints!=0) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">nkpoints [<a href=\"" << url_WEB << "/?nkpoints\">nkpoints</a>]:</h5><span class=\"value\">" << aentry.nkpoints << " </span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.nkpoints_irreducible!=0) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">nkpoints_irreducible [<a href=\"" << url_WEB << "/?nkpoints_irreducible\">nkpoints_irreducible</a>]:</h5><span class=\"value\">" << aentry.nkpoints_irreducible << " </span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.kppra!=0) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">kppra [<a href=\"" << url_WEB << "/?kppra\">kppra</a>]:</h5><span class=\"value\">" << aentry.kppra << " </span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.kpoints.empty()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">kpoints [<a href=\"" << url_WEB << "/?kpoints\">kpoints</a>]:</h5><span class=\"value\">" << aentry.kpoints << " </span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //  oss.precision(3);
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl;  //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Calculation properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::EDATA_ORIG")&& !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << line_rule << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- Warning: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<span class=\"title\"><FONT SIZE=+3 color=red> Warning!</font></span>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<span class=\"title\"><FONT SIZE=+3 color=red> Original and relaxed structures</font></span>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<span class=\"title\"><FONT SIZE=+3 color=red> have different symmetries: listing both.</font></span>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- Warning: END -->" << endl;
//[SC20200327 - OBSOLETE]        oss << line_rule << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // EDATA ORIG/RELAX
//[SC20200327 - OBSOLETE]      if((vflags.flag("FLAG::EDATA_ORIG") || vflags.flag("FLAG::EDATA_RELAX")) && !directory.empty()) {
//[SC20200327 - OBSOLETE]        for(uint i=0;i<=1;i++)  {
//[SC20200327 - OBSOLETE]          if((vflags.flag("FLAG::EDATA_ORIG") && i==0) || (vflags.flag("FLAG::EDATA_RELAX") && i==1)) {
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<!-- EDATA: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            if(i==0) oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Original Structure</h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            if(i==1) oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Relaxed Structure</h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<ht /> " << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class = \"real_space\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            if(i==0) oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Real Space Lattice</h4></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            if(i==1) oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Real Space Lattice</h4></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            vector<string> vline_edata;
//[SC20200327 - OBSOLETE]            if(i==0 && aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT)) aurostd::file2vectorstring(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT,vline_edata);
//[SC20200327 - OBSOLETE]            if(i==1 && aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT)) aurostd::file2vectorstring(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT,vline_edata);
//[SC20200327 - OBSOLETE]            // 
//[SC20200327 - OBSOLETE]            vector<double> abcR(6);
//[SC20200327 - OBSOLETE]            double volumeR,density=0.0,coveraR;
//[SC20200327 - OBSOLETE]            string Crystal_Real_space_Bravais_Lattice_Primitive="",Crystal_Real_space_Lattice_Variation="",Crystal_Real_space_Lattice_System="";
//[SC20200327 - OBSOLETE]            string Crystal_Real_space_Pearson_Symbol="",Crystal_Real_space_Crystal_Family="",Crystal_Real_space_Crystal_System="";
//[SC20200327 - OBSOLETE]            string Crystal_Real_space_Crystal_Class="",Crystal_Real_space_Point_Group_Hermann_Mauguin="",Crystal_Real_space_Point_Group_Schoenflies="";
//[SC20200327 - OBSOLETE]            string Crystal_Real_space_Point_Group_Orbifold="",Crystal_Real_space_Point_Group_Type="",Crystal_Real_space_Point_Group_Order="";
//[SC20200327 - OBSOLETE]            string Crystal_Real_space_Point_Group_Structure="";
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            string Lattice_Real_space_Bravais_Lattice_Primitive="",Lattice_Real_space_Lattice_Variation="",Lattice_Real_space_Lattice_System="";
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            string Superattice_Real_space_Bravais_Superlattice_Primitive="",Superattice_Real_space_Superlattice_Variation="",Superattice_Real_space_Superlattice_System="";
//[SC20200327 - OBSOLETE]            string Superattice_Real_space_Pearson_Symbol_Superlattice="",Reciprocal_lattice_primitive="",Reciprocal_lattice_variation="";
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            vector<double> abcK(6);
//[SC20200327 - OBSOLETE]            double volumeK;
//[SC20200327 - OBSOLETE]            if(vline_edata.size()>0) {
//[SC20200327 - OBSOLETE]              for(uint iline=0;iline<vline_edata.size();iline++) {
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"Real space a b c alpha beta gamma"))
//[SC20200327 - OBSOLETE]                  if(!aurostd::substring2bool(vline_edata.at(iline),"Bohrs/Degs")) {
//[SC20200327 - OBSOLETE]                    aurostd::string2tokens(vline_edata.at(iline),tokens);
//[SC20200327 - OBSOLETE]                    for(uint i=0;i<6;i++) abcR.at(i)=aurostd::string2utype<double>(tokens.at(tokens.size()-6+i));}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"Real space Volume")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline),tokens);
//[SC20200327 - OBSOLETE]                  volumeR=aurostd::string2utype<double>(tokens.at(tokens.size()-1));}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"Real space c/a")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline),tokens);
//[SC20200327 - OBSOLETE]                  coveraR=aurostd::string2utype<double>(tokens.at(tokens.size()-1));}
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"BRAVAIS LATTICE OF THE CRYSTAL")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Crystal_Real_space_Bravais_Lattice_Primitive=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Crystal_Real_space_Lattice_Variation=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Crystal_Real_space_Lattice_System=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+4),tokens,"="); Crystal_Real_space_Pearson_Symbol=tokens.at(tokens.size()-1);}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"POINT GROUP CRYSTAL")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Crystal_Real_space_Crystal_Family=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Crystal_Real_space_Crystal_System=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Crystal_Real_space_Crystal_Class=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+4),tokens,"="); Crystal_Real_space_Point_Group_Hermann_Mauguin=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+5),tokens,"="); Crystal_Real_space_Point_Group_Schoenflies=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+6),tokens,"="); Crystal_Real_space_Point_Group_Orbifold=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+7),tokens,"="); Crystal_Real_space_Point_Group_Type=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+8),tokens,"="); Crystal_Real_space_Point_Group_Order=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+9),tokens,"="); Crystal_Real_space_Point_Group_Structure=tokens.at(tokens.size()-1);}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"BRAVAIS LATTICE OF THE LATTICE")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Lattice_Real_space_Bravais_Lattice_Primitive=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Lattice_Real_space_Lattice_Variation=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Lattice_Real_space_Lattice_System=tokens.at(tokens.size()-1);}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"SUPERLATTICE")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Superattice_Real_space_Bravais_Superlattice_Primitive=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Superattice_Real_space_Superlattice_Variation=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Superattice_Real_space_Superlattice_System=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+4),tokens,"="); Superattice_Real_space_Pearson_Symbol_Superlattice=tokens.at(tokens.size()-1);}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"Reciprocal space a b c alpha beta gamma")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline),tokens);
//[SC20200327 - OBSOLETE]                  for(uint i=0;i<6;i++) abcK.at(i)=aurostd::string2utype<double>(tokens.at(tokens.size()-6+i));}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"Reciprocal space Volume")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline),tokens);
//[SC20200327 - OBSOLETE]                  volumeK=aurostd::string2utype<double>(tokens.at(tokens.size()-1));}
//[SC20200327 - OBSOLETE]                if(aurostd::substring2bool(vline_edata.at(iline),"RECIPROCAL LATTICE")) {
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+7),tokens,"="); Reciprocal_lattice_primitive=tokens.at(tokens.size()-1);
//[SC20200327 - OBSOLETE]                  aurostd::string2tokens(vline_edata.at(iline+8),tokens,"="); Reciprocal_lattice_variation=tokens.at(tokens.size()-1);}
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if(aentry.vcomposition.size()==aentry.vspecies.size()) {
//[SC20200327 - OBSOLETE]              for(uint i=0;i<aentry.vspecies.size();i++) 
//[SC20200327 - OBSOLETE]                density=density+aentry.vcomposition.at(i)*GetAtomMass(aentry.vspecies.at(i))/volumeR*1000.0*1e8*1e8*1e8; //grams/cm^3
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            // Print out structural data    
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<li>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Lattice:</h5>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class=\"lattice_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<table class=\"lattice\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"value\">" << endl << "a=" << abcR.at(0) << "&Aring;" << "&nbsp;" << endl << " b=" << abcR.at(1) << "&Aring;" << "&nbsp;" <<  endl << " c=" << abcR.at(2) << "&Aring;" << "&nbsp;" << endl << " c/a=" << coveraR <<  endl << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"value\">" << endl << "&alpha;=" << abcR.at(3) << "&deg" << endl << " &beta;=" << abcR.at(4) << "&deg" << endl << " &gamma;=" << abcR.at(5) << "&deg" << endl << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - START
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 vector<double> lattice_params(6);
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 vector<string> rtokens;
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 aurostd::string2tokens(aentry.geometry,rtokens,";");
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 for(uint t=0;t<rtokens.size();t++){ lattice_params[t] = aurostd::string2utype<double>(rtokens[t]); }
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 double covera = lattice_params.at(0)/lattice_params.at(2);
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"value\">" << endl << "a=" << lattice_params.at(0) << "&Aring;" << "&nbsp;" << endl << " b=" << lattice_params.at(1) << "&Aring;" << "&nbsp;" <<  endl << " c=" << lattice_params.at(2) << "&Aring;" << "&nbsp;" << endl << " c/a=" << covera <<  endl << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"value\">" << endl << "&alpha;=" << lattice_params.at(3) << "&deg" << endl << " &beta;=" << lattice_params.at(4) << "&deg" << endl << " &gamma;=" << lattice_params.at(5) << "&deg" << endl << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - END
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</table>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</li>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "</div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << volumeR << "&Aring;<sup>3</sup></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - START
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << aentry.volume_cell << "&Aring;<sup>3</sup></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - END
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Number of Atoms per Cell:</h5><span class=\"value\">" << aentry.natoms << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //  if(html && aentry.density==0.0) oss << "<li><span class=\"description\"> Density(calc):</span> " << density << " g/cm<sup>3</sup></li>" << endl;
//[SC20200327 - OBSOLETE]            // if(html && aentry.density>0.0) oss << "<li><span class=\"description\"> Density(entry):</span> " << aentry.density << " g/cm<sup>3</sup></li>" << endl;
//[SC20200327 - OBSOLETE]            if(aentry.density<0.1)  aentry.density=density;
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Density(calc):</h5><span class=\"value\"> " << density << " g/cm<sup>3</sup></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Density(aflowlib):</h5><span class=\"value\"> " << aentry.density << " g/cm<sup>3</sup></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            //	  if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax")) {
//[SC20200327 - OBSOLETE]            //    oss << "<li><span class=\"description\"> Relaxed position (aflowlib/VASP):</span>" << "[<a href=\"" << url_WEB << "/CONTCAR.relax\">POSCAR</a>]" << "</li>" << endl;
//[SC20200327 - OBSOLETE]            // }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.vasp")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/VASP):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/CONTCAR.relax.vasp\"" << html_TAB << ">VASP-POSCAR</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.qe")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/QE):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/CONTCAR.relax.qe\"" << html_TAB << ">QE-GEOMETRY</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.abinit")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/ABINIT):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/CONTCAR.relax.abinit\"" << html_TAB << ">ABINIT-GEOMETRY</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.aims")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/AIMS):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/CONTCAR.relax.aims\"" << html_TAB << ">AIMS-GEOMETRY</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"INCAR.relax")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> INCAR for relax calculation (aflowlib/VASP):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/INCAR.relax\"" << html_TAB << ">INCAR.relax</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if(aurostd::substring2bool(aentry.vfiles_WEB,"INCAR.static")) {
//[SC20200327 - OBSOLETE]              oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> INCAR for static calculation (aflowlib/VASP):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                << "[<a href=\"" << url_WEB << "/INCAR.static\"" << html_TAB << ">INCAR.static</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"INCAR.bands")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> INCAR for bands calculation (aflowlib/VASP):</h5><span class=\"value\">"
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/INCAR.bands\"" << html_TAB << ">INCAR.bands</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"KPOINTS.relax")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> KPOINTS for relax calculation (aflowlib/VASP):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/KPOINTS.relax\"" << html_TAB << ">KPOINTS.relax</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"KPOINTS.static")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> KPOINTS for static calculation (aflowlib/VASP):</h5><span class=\"value\">"
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/KPOINTS.static\"" << html_TAB << ">KPOINTS.static</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,"KPOINTS.bands")) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> KPOINTS for bands calculation (aflowlib/VASP):</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/KPOINTS.bands\"" << html_TAB << ">KPOINTS.bands</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]            if(aurostd::substring2bool(aentry.vfiles_WEB,DEFAULT_FILE_EDATA_ORIG_OUT)) {
//[SC20200327 - OBSOLETE]              oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Extended crystallographic data for original structure:</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                << "[<a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_ORIG_OUT << "\"" << html_TAB << ">" << DEFAULT_FILE_EDATA_ORIG_OUT << "</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,DEFAULT_FILE_EDATA_RELAX_OUT)) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Extended crystallographic data for relaxed structure:</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_RELAX_OUT << "\"" << html_TAB << ">" << DEFAULT_FILE_EDATA_RELAX_OUT << "</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            if (i==1){  //GG20170714 - removed relaxed data from orig
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB,DEFAULT_FILE_EDATA_BANDS_OUT)) {
//[SC20200327 - OBSOLETE]                oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Extended crystallographic data for band-structure:</h5><span class=\"value\">" 
//[SC20200327 - OBSOLETE]                  << "[<a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_BANDS_OUT << "\"" << html_TAB << ">" << DEFAULT_FILE_EDATA_BANDS_OUT << "</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              }
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class=\"space_group\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Bravais Lattice of the Crystal" << aflow_sym_readme << art135_link << "</h4></div>" << endl; //PC20180620 //JPO20180731  //CO20180817
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Space Group Number:</h5><span class=\"value\">" << aentry.sg2 << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Pearson Symbol:</h5><span class=\"value\">" << Crystal_Real_space_Pearson_Symbol << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive:</h5><span class=\"value\">" << Crystal_Real_space_Bravais_Lattice_Primitive << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << Crystal_Real_space_Lattice_Variation << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << Crystal_Real_space_Lattice_System << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class=\"point_group\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Point Group of the Crystal" << aflow_sym_readme << art135_link << "</h4></div>" << endl; //PC20180620 //JPO20180731 //CO20180817
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Family:</h5><span class=\"value\">" << Crystal_Real_space_Crystal_Family << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal System:</h5><span class=\"value\">" << Crystal_Real_space_Crystal_System << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Class:</h5><span class=\"value\">" << Crystal_Real_space_Crystal_Class << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //	oss << "<li><span class=\"description\"> Point Group (Hermann Mauguin):</span>" << Crystal_Real_space_Point_Group_Hermann_Mauguin << "</li><!br>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Herm. Maug.):</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Hermann_Mauguin << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Schoenflies):</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Schoenflies << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Orbifold:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Orbifold << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180827 - new point group type output - START
//[SC20200327 - OBSOLETE]            string point_group_type = Crystal_Real_space_Point_Group_Type;
//[SC20200327 - OBSOLETE]            if(aurostd::RemoveWhiteSpaces(point_group_type) == "-" || aurostd::RemoveWhiteSpaces(point_group_type) == "none"){
//[SC20200327 - OBSOLETE]              point_group_type = "non-centrosymmetric, non-enantiomorphic, non-polar";
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << point_group_type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180827 - new point group type output - END
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180827 [OBSOLETE] oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Order:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Order << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Structure:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Structure << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - START
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Family:</h5><span class=\"value\">" << aentry.crystal_family << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal System:</h5><span class=\"value\">" << aentry.crystal_system << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Class:</h5><span class=\"value\">" << aentry.crystal_class << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Herm. Maug.):</h5><span class=\"value\">" << aentry.point_group_Hermann_Mauguin << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Schoenflies):</h5><span class=\"value\">" << aentry.point_group_Schoenflies << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Orbifold:</h5><span class=\"value\">" << aentry.point_group_orbifold << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << aentry.point_group_type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Order:</h5><span class=\"value\">" << aentry.point_group_order << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Structure:</h5><span class=\"value\">" << aentry.point_group_structure << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - END
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class=\"bravais_lattice\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Bravais Lattice of the Lattice" << aflow_sym_readme << art135_link << "</h4></div>" << endl; //PC20180620 //JPO20180731 //CO20180817
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive</h5><span class=\"value\">" << Lattice_Real_space_Bravais_Lattice_Primitive << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_Variation << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_System << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - START
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive</h5><span class=\"value\">" << aentry.Bravais_lattice_lattice_type << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << aentry.Bravais_lattice_lattice_variation_type << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << aentry.Bravais_lattice_lattice_system << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - END
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            if(aurostd::substring2bool(aentry.vfiles_WEB,label+"_BZ.png")) {
//[SC20200327 - OBSOLETE]              oss << "<div class=\"container__cell\"><div class=\"container__card--img\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              oss << "<h5 class=\"value-name\"> Brillouin Zone " << art058_link<< "</h5>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              // [OBSOLETE] oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              // [OBSOLETE] oss << "<div class=\"picture_BZ\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]              // [OBSOLETE] oss << "<img class=\"pic_BZ\" src=\"../SCIENCE/images/brillouin/" << aurostd::RemoveWhiteSpaces(Lattice_Real_space_Lattice_Variation) << ".PNG\" alt=\"Brillouin Zone of " << label << "\" />" << endl; //CO20170621 - relative path
//[SC20200327 - OBSOLETE]              oss << "<img class=\"BZ-img\" src=\"http://aflowlib.duke.edu/SCIENCE/images/brillouin/" << aurostd::RemoveWhiteSpaces(Lattice_Real_space_Lattice_Variation) << ".PNG\" alt=\"Brillouin Zone of " << label << "\" />" << endl; //CO20170621 - abs path WORKS  //JPO20180731
//[SC20200327 - OBSOLETE]              oss << "</div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]// [OBSOLETE] oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]// [OBSOLETE] oss << "<div class=\"bravais_lattice\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Bravais Lattice of the Lattice" << aflow_sym_readme << art135_link << "</h4></div>" << endl; //PC20180620 //JPO20180731 //CO20180817
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]// [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive</h5><span class=\"value\">" << Lattice_Real_space_Bravais_Lattice_Primitive << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_Variation << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_System << "</span></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]// [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            //[MOVED UP JPO20180731]oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class=\"superlattice\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Superlattice" << aflow_sym_readme << art135_link << "</h4></div>" << endl; //PC20180620 //JPO20180731 //CO20180817
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Primitive unit cell:</h5><span class=\"value\">" << Superattice_Real_space_Bravais_Superlattice_Primitive << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Variation:</h5><span class=\"value\">" << Superattice_Real_space_Superlattice_Variation << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Lattice System :</h5><span class=\"value\">" << Superattice_Real_space_Superlattice_System << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Pearson Symbol:</h5><span class=\"value\">" << Superattice_Real_space_Pearson_Symbol_Superlattice << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - START
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Primitive unit cell:</h5><span class=\"value\">" << aentry.Bravais_superlattice_lattice_type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Variation:</h5><span class=\"value\">" << aentry.Bravais_superlattice_lattice_variation_type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Lattice System :</h5><span class=\"value\">" << aentry.Bravais_superlattice_lattice_system << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Pearson Symbol:</h5><span class=\"value\">" << aentry.Pearson_symbol_superlattice << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - END
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class=\"reciprocal\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Reciprocal Space Lattice </h4></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Reciprocal Lattices:</h5>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<div class=\"lattice_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<table class=\"reciprocal_lattice\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"value\">" << endl << "a=" << abcK.at(0) << "&Aring;<sup>-1</sup>" << endl << "b=" << abcK.at(1) << "&Aring;<sup>-1</sup>" << endl << "c=" << abcK.at(2) << "&Aring;<sup>-1</sup>" << endl << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"value\">" << "&alpha;=" << abcK.at(3) << "&deg;" << endl << "&beta;=" << abcK.at(4) << "&deg;" << endl << "&gamma;=" << abcK.at(5) << "&deg;" << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - START
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 vector<double> reciprocal_lattice_params(6);
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 vector<string> ktokens;
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 aurostd::string2tokens(aentry.reciprocal_geometry,ktokens,";");
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 for(uint t=0;t<ktokens.size();t++){ reciprocal_lattice_params[t] = aurostd::string2utype<double>(ktokens[t]); }
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"value\">" << endl << "a=" << reciprocal_lattice_params.at(0) << "&Aring;<sup>-1</sup>" << endl << "b=" << reciprocal_lattice_params.at(1) << "&Aring;<sup>-1</sup>" << endl << "c=" << reciprocal_lattice_params.at(2) << "&Aring;<sup>-1</sup>" << endl << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"value\">" << "&alpha;=" << reciprocal_lattice_params.at(3) << "&deg;" << endl << "&beta;=" << reciprocal_lattice_params.at(4) << "&deg;" << endl << "&gamma;=" << reciprocal_lattice_params.at(5) << "&deg;" << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - END
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "<tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</table>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "</div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << volumeK << " &Aring;<sup>-3</sup></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Lattice Primitive:</h5><span class=\"value\">" << Reciprocal_lattice_primitive << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<div class=\"container__cell\"><div class =\"container__card\"><h5 class=\"value-name\"> Lattice Variation:</h5><span class=\"value\">" << Reciprocal_lattice_variation << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - START
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << aentry.reciprocal_volume_cell << " &Aring;<sup>-3</sup></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Lattice Primitive:</h5><span class=\"value\">" << aentry.reciprocal_lattice_type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 oss << "<div class=\"container__cell\"><div class =\"container__card\"><h5 class=\"value-name\"> Lattice Variation:</h5><span class=\"value\">" << aentry.reciprocal_lattice_variation_type << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            // [OBSOLETE] DX20180824 - will replace lines above when the database is sufficiently populated - END
//[SC20200327 - OBSOLETE]            // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "</div>" << endl;  //JPO20180731
//[SC20200327 - OBSOLETE]            oss << "<!-- EDATA: END -->" << endl;
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // THERMODYNAMICS 
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::THERMODYNAMICS") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- Thermodynamics properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"Thermodynamics\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\">Thermodynamics Properties</h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //   oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //     oss << "<br>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<span class=\"Thermodynamics_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Formation Enthalpy_cell:</h5><span class=\"value\">" << (aentry.enthalpy_formation_cell!=0.0?aurostd::utype2string(aentry.enthalpy_formation_cell,5)+" (eV) ":"unavailable") << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Formation Enthalpy_atom:</h5><span class=\"value\">" << (aentry.enthalpy_formation_atom!=0.0?aurostd::utype2string(aentry.enthalpy_formation_atom,5)+" (eV/atom) ":"unavailable") << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"><i>Ab-initio</i> energy_cell:</h5><span class=\"value\">" << (aentry.energy_cell!=0.0?aurostd::utype2string(aentry.energy_cell,5)+" (eV) ":"unavailable") << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"><i>Ab-initio</i> energy_atom:</h5><span class=\"value\">" << (aentry.energy_cell!=0.0?aurostd::utype2string(aentry.energy_atom,5)+" (eV/atom) ":"unavailable") << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Thermodynamics properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // AGL
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::AGL") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- AGL properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"AGL\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\">AGL Properties (Aflow Gibbs Library)" << aflow_agl_readme << art115_link << "</h1></div>" << endl; //PC20180620 //JPO20180731  //CO20180817
//[SC20200327 - OBSOLETE]        //   oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //     oss << "<br>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<span class=\"AGL_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.agl.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Output" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/aflow.agl.out\"" << html_TAB << ">aflow.agl.out</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Complete Output" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AGL.out\"" << html_TAB << ">AGL.out</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_energies_temperature.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Energy versus Temperature" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AGL_energies_temperature.out\"" << html_TAB << ">AGL_energies_temperature.out</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_thermal_properties_temperature.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Thermal Properties versus Temperature" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AGL_thermal_properties_temperature.out\"" << html_TAB << ">AGL_thermal_properties_temperature.out</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_Hugoniot.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Hugoniot Relation" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AGL_Hugoniot.out\"" << html_TAB << ">AGL_Hugoniot.out</a>]" << "</span></div></div>" << endl; //XCT 20181212
//[SC20200327 - OBSOLETE]        if(aentry.agl_thermal_conductivity_300K<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_thermal_conductivity_300K_td\">AGL Thermal Conductivity at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_thermal_conductivity_300K << " (W/m*K)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_debye<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_debye_td\">AGL Debye Temperature" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_debye << " (K)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_acoustic_debye<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_acoustic_debye_td\">AGL Debye Acoustic Temperature" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_acoustic_debye << " (K)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_gruneisen<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_gruneisen_td\">AGL Gruneisen parameter" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_gruneisen << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_heat_capacity_Cv_300K<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_heat_capacity_Cv_300K_td\">AGL Specific Heat Cv at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_heat_capacity_Cv_300K << " (kB/cell)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_heat_capacity_Cp_300K<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_heat_capacity_Cp_300K_td\">AGL Specific Heat Cp at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_heat_capacity_Cp_300K << " (kB/cell)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_thermal_expansion_300K<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_thermal_expansion_300K_td\">AGL Thermal Expansion at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_thermal_expansion_300K << " (1/K)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_bulk_modulus_static_300K<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_bulk_modulus_static_300K_td\">AGL Bulk Modulus Static at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_bulk_modulus_static_300K << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_bulk_modulus_isothermal_300K<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_bulk_modulus_isothermal_300K_td\">AGL Bulk Modulus Isothermal at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_bulk_modulus_isothermal_300K << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.agl_poisson_ratio_source.size()) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_poisson_ratio_source_td\">AGL Poisson Ratio Source" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_poisson_ratio_source << " </span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.agl_vibrational_free_energy_300K_cell<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_vibrational_free_energy_300K_cell_td\">AGL Vibrational Free Energy per cell at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_vibrational_free_energy_300K_cell << " (meV/cell)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.agl_vibrational_free_energy_300K_atom<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_vibrational_free_energy_300K_atom_td\">AGL Vibrational Free Energy per atom at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_vibrational_free_energy_300K_atom << " (meV/atom)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.agl_vibrational_entropy_300K_cell<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_vibrational_entropy_300K_cell_td\">AGL Vibrational Entropy per cell at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_vibrational_entropy_300K_cell << " (meV/cell*K)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.agl_vibrational_entropy_300K_atom<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_vibrational_entropy_300K_atom_td\">AGL Vibrational Entropy per atom at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_vibrational_entropy_300K_atom << " (meV/atom*K)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- AGL properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // AEL
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::AEL") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- AEL properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"AEL\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\">AEL Properties (Aflow Elastic Library)</font>" << aflow_ael_readme << art096_link << "</h1></div>" << endl; //PC20180620 //JPO20180731  //CO20180817
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<span class=\"AEL_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.ael.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AEL Output" << art100_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/aflow.ael.out\"" << html_TAB << ">aflow.ael.out</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Elastic_constants.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AEL Elastic Constants (stiffness tensor)" << art100_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AEL_Elastic_constants.out\"" << html_TAB << ">AEL_Elastic_constants.out</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Compliance_tensor.out")) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AEL Compliance Constants (compliance tensor)" << art100_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AEL_Compliance_tensor.out\"" << html_TAB << ">AEL_Compliance_tensor.out</a>]" << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_poisson_ratio<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_poisson_ratio_td\">AEL Poisson Ratio" << art100_link << ":</h5><span class=\"value\"> " << aentry.ael_poisson_ratio << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_bulk_modulus_voigt<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_bulk_modulus_voigt_td\">AEL Bulk Modulus Voigt" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_bulk_modulus_voigt << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_bulk_modulus_reuss<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_bulk_modulus_reuss_td\">AEL Bulk Modulus Reuss" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_bulk_modulus_reuss << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_bulk_modulus_vrh<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_bulk_modulus_vrh_td\">AEL Bulk Modulus VRH" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_bulk_modulus_vrh << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_shear_modulus_voigt<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_shear_modulus_voigt_td\">AEL Shear Modulus Voigt" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_shear_modulus_voigt << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_shear_modulus_reuss<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_shear_modulus_reuss_td\">AEL Shear Modulus Reuss" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_shear_modulus_reuss << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_shear_modulus_vrh<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_shear_modulus_vrh_td\">AEL Shear Modulus VRH" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_shear_modulus_vrh << " (GPa)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.ael_elastic_anisotropy<AUROSTD_NAN) //CO20181129
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_elastic_anisotropy_td\">AEL Elastic Anisotropy" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_elastic_anisotropy << "</span></div></div>" << endl; //JPO20180731 //CO20181129
//[SC20200327 - OBSOLETE]        if(aentry.ael_youngs_modulus_vrh<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_youngs_modulus_vrh_td\">AEL Young's Modulus VRH" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_youngs_modulus_vrh << " (GPa)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.ael_speed_sound_transverse<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_speed_sound_transverse_td\">AEL Speed of Sound in Transverse Direction" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_speed_sound_transverse << " (m/s)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.ael_speed_sound_longitudinal<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_speed_sound_longitudinal_td\">AEL Speed of Sound in Longitudinal Direction" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_speed_sound_longitudinal << " (m/s)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.ael_speed_sound_average<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_speed_sound_average_td\">AEL Speed of Sound, Average" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_speed_sound_average << " (m/s)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.ael_pughs_modulus_ratio<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_pughs_modulus_ratio_td\">AEL Pugh's Modulus Ratio" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_pughs_modulus_ratio << " </span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.ael_debye_temperature<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_debye_temperature_td\">AEL Debye Temperature" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_debye_temperature << " (K)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.ael_applied_pressure<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_applied_pressure_td\">AEL Applied Pressure" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_applied_pressure << " (GPa)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        if(aentry.ael_average_external_pressure<AUROSTD_NAN) 
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_average_external_pressure_td\">AEL Average External Pressure" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_average_external_pressure << " (GPa)</span></div></div>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- AEL properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // BADER
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::BADER") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- Bader properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"Bader\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\">Bader Atoms in Molecules Properties</font>" << "</h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<span class=\"title\"><FONT SIZE=+3>Bader Atoms in Molecules Properties</font>" << "<a href=http://materials.duke.edu/AFLOW/README_AFLOW_BADER.TXT>[info]</a>" << "</span>" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<span class=\"Bader_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        string abader_out=aentry.prototype+"_abader.out";
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.vfiles_WEB,abader_out)) {oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bader Output" << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/" << abader_out << "\"" << html_TAB << ">" << abader_out << "</a>]" << "</span></div></div>" << endl;} //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<li><span class=\"description\"> Bader Output" << art100_link << ":</span>" << "[<a href=\"" << url_WEB << "/" << abader_out << "\"" << html_TAB << ">" << abader_out << "</a>]" << "</li><!br>" << endl;
//[SC20200327 - OBSOLETE]        string bader_net_charges_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bader_net_charges target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]        string bader_atomic_volumes_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bader_atomic_volumes target=\"_blank\"><font color=black><i>info</i></font></a>]";
//[SC20200327 - OBSOLETE]        if(!aentry.vbader_net_charges.empty()) {
//[SC20200327 - OBSOLETE]          aurostd::string2tokens(aentry.species,tokens,",");  //tokens=species
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"bader_net_charges_td\">Net Charges" << bader_net_charges_wiki_link << ":</h5><span class=\"value\"> "; //JPO20180731
//[SC20200327 - OBSOLETE]          atomCOUNT=0;
//[SC20200327 - OBSOLETE]          for(uint i=0;i<tokens.size();i++) {
//[SC20200327 - OBSOLETE]            oss << tokens.at(i) << "={";
//[SC20200327 - OBSOLETE]            for(uint j=0;j<(uint)aentry.vcomposition.at(i);j++) {
//[SC20200327 - OBSOLETE]              num_prec.str("");
//[SC20200327 - OBSOLETE]              num_prec << std::fixed << setprecision(4) << aentry.vbader_net_charges.at(atomCOUNT++);
//[SC20200327 - OBSOLETE]              oss << num_prec.str();
//[SC20200327 - OBSOLETE]              if(j<(uint)aentry.vcomposition.at(i)-1) {oss << ", ";}
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]            if(i<tokens.size()-1) {oss << "}; ";} else {oss << "}";}
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          oss << " (electrons)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        if(!aentry.bader_atomic_volumes.empty()) {
//[SC20200327 - OBSOLETE]          aurostd::string2tokens(aentry.species,tokens,",");  //tokens=species
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"bader_atomic_volumes_td\">Atomic Volumes" << bader_atomic_volumes_wiki_link << ":</h5><span class=\"value\">"; //JPO20180731
//[SC20200327 - OBSOLETE]          atomCOUNT=0;
//[SC20200327 - OBSOLETE]          for(uint i=0;i<tokens.size();i++) {
//[SC20200327 - OBSOLETE]            oss << tokens.at(i) << "={";
//[SC20200327 - OBSOLETE]            for(uint j=0;j<(uint)aentry.vcomposition.at(i);j++) {
//[SC20200327 - OBSOLETE]              num_prec.str("");
//[SC20200327 - OBSOLETE]              num_prec << std::fixed << setprecision(4) << aentry.vbader_atomic_volumes.at(atomCOUNT++);
//[SC20200327 - OBSOLETE]              oss << num_prec.str();
//[SC20200327 - OBSOLETE]              if(j<(uint)aentry.vcomposition.at(i)-1) {oss << ", ";}
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]            if(i<tokens.size()-1) {oss << "}; ";} else {oss << "}";}
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]          oss << " (Angst<sup>3</sup>)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Bader properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // MAGNETIC
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::MAGNETIC") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- Magnetic properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"Magnetic\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\">Magnetic Properties </h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<!span class=\"Magnetic_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<li><span class=\"description\" id=\"spin_cell_td\">Magnetic Moment:</span> " << aentry.spin_cell << " &mu;<sub>B</sub></li><!br>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<li><span class=\"description\" id=\"spin_atom_td\">Magnetic Moment/atom:</span> " << aentry.spin_atom << " &mu;<sub>B</sub>/atom</li><!br>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aurostd::substring2bool(aentry.loop,"magnetic")) {
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spinD_td\">Spin Decomposition per atoms:</h5><span class=\"value\"> {" << aentry.spinD << "} &mu;<sub>B</sub></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spinF_td\">Spin Polarization (E<sub>F</sub>):</h5><span class=\"value\">" << aentry.spinF << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spin_cell_td\">Magnetic Moment:</h5><span class=\"value\">" << aentry.spin_cell << " &mu;<sub>B</sub></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spin_atom_td\">Magnetic Moment/atom:</h5><span class=\"value\">" << aentry.spin_atom << " &mu;<sub>B</sub>/atom</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</ul>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Magnetic properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // SCINTILLATION
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::SCINTILLATION") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        double scintillation_attenuation_length=GetCompoundAttenuationLength(aentry.vspecies,aentry.vcomposition,aentry.density);
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Scintillation properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"scintillation\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //    oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\">Scintillation Properties</h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        //   oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<br>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<table class=\"scintillation_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.scintillation_attenuation_length>0.0) {
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\">" << endl << "<h5 class=\"value-name\">Attenuation Length" << art064_link << ":</h5><span class=\"value\">" << aentry.scintillation_attenuation_length << " cm" << endl << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        } else {
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\">" << endl << "<h5 class=\"value-name\">Attenuation Length" << art064_link << ":</h5><span class=\"value\">" << scintillation_attenuation_length << " cm" << endl << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</table>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Scintillation properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // ELECTRONIC BANDS
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ELECTRONIC") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << line_rule << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<!-- Electronic properties: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        //   oss << "<hr />" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"electronic\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Electronic Properties </h1></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<table class=\"electronic_table\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\">Spin Decomposition per atoms:</h5><span class=\"value\"> {" << aentry.spinD << "} &mu;<sub>B</sub> </span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Spin Polarization (E<sub>F</sub>):</h5><span class=\"value\"> " << aentry.spinF << " </span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"band_gap_td\"><h5 class=\"value-name\">Band Gap:</h5><span class=\"value\">" << aentry.Egap << " eV "
//[SC20200327 - OBSOLETE]          << (aentry.Egap_type=="insulator_direct"?"(insulator)":"") << (aentry.Egap_type=="insulator_indirect"?"(insulator)":"") << (aentry.Egap_type=="metal"?"(metal)":"") << "</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"fit_band_gap_td\"><h5 class=\"value-name\">Fit Band Gap:</h5><span class=\"value\">" << aentry.Egap_fit << " eV</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</tr>" << endl;  //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_m_atom_td\"><h5 class=\"value-name\">Magnetic Moment:</h5><span class=\"value\">" << aentry.spin_cell << " &mu;<sub>B</sub></span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_m_atom_td\"><h5 class=\"value-name\">Magnetic Moment/atom:</h5><span class=\"value\">" << aentry.spin_atom << " &mu;<sub>B</sub>/atom</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // (No values available - do not display yet) oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_e_mass_td\"><h5 class=\"value-name\">Electron Mass(FIX):</h5><span class=\"value\"> XXX (m<sub>0</sub>)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // (No values available - do not display yet) oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_hole_mass_td\"><h5 class=\"value-name\">Hole Mass(FIX):</h5><span class=\"value\"> XXX (m<sub>0</sub>)</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.Egap_type=="insulator_direct")
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_band_gap_type_td\"><h5 class=\"value-name\">Band Gap Type:</h5><span class=\"value\">  Direct</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        if(aentry.Egap_type=="insulator_indirect")
//[SC20200327 - OBSOLETE]          oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_band_gap_type_td\"><h5 class=\"value-name\">Band Gap Type:</h5><span class=\"value\">  Indirect</span></div></div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<td class=\"electronic_table_td\"><span class=\"table_description\">Spin Polarization (E<sub>F</sub>):</span> " << aentry.spinF << " </td>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<td class=\"electronic_table_td\"><span class=\"table_description\">Spin Decomposition per atoms:</span> {" << aentry.spinD << "} &mu;<sub>B</sub> </td>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</tr>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</tbody>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "</table>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<ul>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<li>" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"container-subtitle\">" << endl;
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<h4 class=\"section-subtitle\">Band Structure</h4></div>" << endl;
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        // ****************************************************************************
//[SC20200327 - OBSOLETE]        // INTERACTIVE BANDS PLOT
//[SC20200327 - OBSOLETE]        //GG
//[SC20200327 - OBSOLETE]        // [OBSOLETE] oss << "<div class=\"container__cell--full\" >" << endl;  //PC20180515 //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"flex-container\">" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<span class=\"pic_description_band\">Band Structure:</span>" << endl; //JPO20180731
//[SC20200327 - OBSOLETE]        oss << "<div class=\"DosOptions\">" 				<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<div> Zoom/Pan type:" << endl;  //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<select id=\"zoomOptions\"> "								<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<option value=\"both\">Both X & Y</option>"				<< endl;
//[SC20200327 - OBSOLETE]        oss << "<option value=\"xOnly\">X Only</option>"				<< endl;
//[SC20200327 - OBSOLETE]        oss << "<option value=\"yOnly\">Y Only</option>"				<< endl;
//[SC20200327 - OBSOLETE]        oss << "</select></div>" << endl; //PC20180515
//[SC20200327 - OBSOLETE]        //if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spinF),0.0,_ZERO_TOL_))) //PC20180515
//[SC20200327 - OBSOLETE]        if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spin_atom),0.0,_ZERO_TOL_))) //PC20180515
//[SC20200327 - OBSOLETE]        { //CO20200106 - patching for auto-indenting
//[SC20200327 - OBSOLETE]          oss << "<div>Majority/Minority Spin Selection:**" << endl;  //PC20180525
//[SC20200327 - OBSOLETE]          oss << "<select id =\"spinBandsOptions\">" << endl; //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<option value=\"bothS\">Both Spins</option>"				<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<option value=\"majority\">Majority Spin</option>"				<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<option value=\"minority\">Minority Spin</option>"				<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "</select></div>"								<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        oss << "<div class=\"reset\">Reset Zoom</div>" << endl; //PC20180515
//[SC20200327 - OBSOLETE]        oss << "</div>"				<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]        if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spin_atom),0.0,_ZERO_TOL_))){  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<div class=\"star\">**For smoother tracing of bands </div>" << endl; //PC20180525
//[SC20200327 - OBSOLETE]        } //PC20180525
//[SC20200327 - OBSOLETE]        oss << "<div class=\"plots\">" << endl; //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<div class=\"Bands_plot\">" << endl;  //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<svg id=\"bands_wrapper\" width=\"900\" height=\"550\"></svg>"						<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]        if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spin_atom),0.0,_ZERO_TOL_))){  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<div id=\"bandLegend\" class=\"legend\">"				<< endl;
//[SC20200327 - OBSOLETE]          oss << "<text class=\"legendText \">Majority Spin"				<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<line class=\"legendLine\" style=\"border-color:black;\"></line>"	        << endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "</text>" 								<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<text class=\"legendText \">Minority Spin" 				<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "<line class=\"legendLine min\"></line>" 					<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]          oss << "</text></div>"							<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]        } //PC20180515
//[SC20200327 - OBSOLETE]        oss << "</div>"  << endl; //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<div class=\"Dos_plot\">" << endl;  //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<svg id=\"dos_wrapper\" width=\"300px\" height=\"550px\" ></svg>" << endl;  //PC20180515
//[SC20200327 - OBSOLETE]        oss << " <div id=\"dosLegend\" class=\"legend\" >" << endl; //PC20180515
//[SC20200327 - OBSOLETE]        oss << "<text id=\"dosText\"></text>"						<< endl;  //PC20180515
//[SC20200327 - OBSOLETE]        oss << "</div>" 								<< endl;
//[SC20200327 - OBSOLETE]        oss << "</div></div></div>" << endl;  //PC20180515
//[SC20200327 - OBSOLETE]        //[OBSOLETE PC20180515]oss << "<div class=\"legendLine min\"></div>" 					<< endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE PC20180515]oss << "</div></div></svg>"							<< endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE PC20180515]oss << "<svg id=\"dos_wrapper\">"						<< endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE PC20180515]oss << "<div id=\"dosText\"></div>"						<< endl;
//[SC20200327 - OBSOLETE]        //[OBSOLETE PC20180515]oss << "<div id=\"dosLegend\" class=\"legend \"></div></svg></div></div>"	<< endl;
//[SC20200327 - OBSOLETE]        // ********************************************************************************** 
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]        if(aentry.vfiles_WEB.size()>0) 
//[SC20200327 - OBSOLETE]          for(uint i=0;i<aentry.vfiles_WEB.size();i++)
//[SC20200327 - OBSOLETE]            if((aentry.vfiles_WEB.at(i)==label+".png") || (aentry.vfiles_WEB.at(i) == label + "_banddos.png"))  //ME20190621 - include new file naming convention 
//[SC20200327 - OBSOLETE]              oss << "<div class = \"picture_band\" id=\"band_dos_pic\"><a id=\"imgPopup\"><img class=\"pic\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(i) << "\" alt=\"Band Structure of " << label << "\" style='display:block;background:white; margin: 0 auto;' /></a></div>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<div id=\"small_figure\" style=\'display:flex;justify-content:center;\'>" << endl;
//[SC20200327 - OBSOLETE]        if(aentry.vfiles_WEB.size()>0) 
//[SC20200327 - OBSOLETE]          for(uint i=0;i<aentry.vfiles_WEB.size();i++)
//[SC20200327 - OBSOLETE]            if((aurostd::substring2bool(aentry.vfiles_WEB.at(i),"_PEDOS_") || aurostd::substring2bool(aentry.vfiles_WEB.at(i), "_dos_")) && aurostd::substring2bool(aentry.vfiles_WEB.at(i),".png"))  //ME20190621 - include new file name convention
//[SC20200327 - OBSOLETE]              oss << "<img class=\"pic_small\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(i) << "\" alt=\"Band Structure of pdos_" << label << "\" style='background:white; margin: 0 auto;height:0%;' >" << endl;
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl;
//[SC20200327 - OBSOLETE]        oss << "</li></ul>" << endl;
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- Electronic properties: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // ELECTRONIC POPUP
//[SC20200327 - OBSOLETE]      if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ELECTRONIC") && !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << "<!-- Electronic popup: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "</div> <!close global>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<div id=\"popupImg\">  <! make popup>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<a id=\"popupImgClose\">" << endl;
//[SC20200327 - OBSOLETE]        oss << "<img id=\"closeButton\" src=\"img/closeButton.png\" alt=\"close\" ></a>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<a id=\"backwardArrow\">" << endl;
//[SC20200327 - OBSOLETE]        oss << "<img src=\"img/backwardArrow.png\" alt=\"backward arrow\" ></a>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<ul class=\"pic_popup_list\">" << endl;
//[SC20200327 - OBSOLETE]        //ME20190621 BEGIN
//[SC20200327 - OBSOLETE]        // Include new file naming convention
//[SC20200327 - OBSOLETE]        if(aentry.vfiles_WEB.size()>0) {
//[SC20200327 - OBSOLETE]          for(uint i=0;i<aentry.vfiles_WEB.size();i++) {
//[SC20200327 - OBSOLETE]            if((aentry.vfiles_WEB.at(i)==label+".png") || (aentry.vfiles_WEB.at(i) == label + "_banddos.png")) {
//[SC20200327 - OBSOLETE]              oss << "<li class=\"pic_popup\"><img class=\"pic_large\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(i) << "\" alt=\"Band Structure of " << label << "\"  ></li>" << endl;
//[SC20200327 - OBSOLETE]              break;
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        if(aentry.vfiles_WEB.size()>0) {
//[SC20200327 - OBSOLETE]          for(uint i=0;i<aentry.vfiles_WEB.size();i++) {
//[SC20200327 - OBSOLETE]            if((aentry.vfiles_WEB.at(i)==label+"_DOS.png") || (aentry.vfiles_WEB.at(i) == label + "_dos.png")) {
//[SC20200327 - OBSOLETE]              oss << "<li class=\"pic_popup\"><img class=\"pic_large\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(i) << "\" alt=\"DOS of " << label << "\"  ></li>" << endl;
//[SC20200327 - OBSOLETE]              break;
//[SC20200327 - OBSOLETE]            }
//[SC20200327 - OBSOLETE]          }
//[SC20200327 - OBSOLETE]        }
//[SC20200327 - OBSOLETE]        //ME20190621 END
//[SC20200327 - OBSOLETE]        if(aentry.vfiles_WEB.size()>0) 
//[SC20200327 - OBSOLETE]          for(uint iline=0;iline<aentry.vfiles_WEB.size();iline++)
//[SC20200327 - OBSOLETE]            if(aurostd::substring2bool(aentry.vfiles_WEB.at(iline),"PEDOS") || aurostd::substring2bool(aentry.vfiles_WEB.at(iline), "_dos_"))  //ME20190621 - include new file name convention
//[SC20200327 - OBSOLETE]              if(aurostd::substring2bool(aentry.vfiles_WEB.at(iline),"png"))
//[SC20200327 - OBSOLETE]                oss << "<li class=\"pic_popup\"><img class=\"pic_large\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(iline) << "\" alt=\"Band Structure of pdos_" << label << "\"  ></li>" << endl;
//[SC20200327 - OBSOLETE]        oss << "</ul>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<a id=\"forwardArrow\"><img src=\"img/forwardArrow.png\" alt=\"forward arrow\" ></a>" << endl;
//[SC20200327 - OBSOLETE]        oss << "</div>" << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- Electronic popup: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      oss << "<!-- Downloadable Files: BEGIN -->" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<div class=\"container\">" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Downloadable Files </h1></div>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<div class=\"container__cell--full\"><div class=\"container__card\">" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<ul class=\"file-list\">" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/\"" << html_TAB << " download=\"aflowlib.out\" >aflowlib.out</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/?format=json\"" << html_TAB << "download=\"aflowlib.json\">aflowlib.json</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.vasp\"" << html_TAB << "download=\"CONTCAR.relax.vasp\">Relaxed Position (aflowlib/VASP)</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.qe\"" << html_TAB << " download=\"CONTCAR.relax.ge\">Relaxed Position (aflowlib/QE)</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.abinit\"" << html_TAB << "download=\"CONTCAR.relax.abinit\">Relaxed Position (aflowlib/AIMS)</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.aims\"" << html_TAB << "download=\"CONTCAR.relax.aims\">Relaxed Position (aflowlib/AIMS)</a>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/INCAR.relax\"" << html_TAB << "download=\"INCAR.relax\">INCAR for relaxed calculation</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/INCAR.static\"" << html_TAB << "download=\"INCAR.static\">INCAR for static calculation</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/INCAR.bands\"" << html_TAB << "download=\"INCAR.bands\">INCAR for bands calculation</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/KPOINTS.relax\"" << html_TAB << "download=\"KPOINTS.relax\">KPOINTS for relaxed calculation</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/KPOINTS.static\"" << html_TAB << "download=\"KPOINTS.static\">KPOINTS for static calculation</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/KPOINTS.bands\"" << html_TAB << "download=\"KPOINTS.bands\">KPOINTS for bands calculation</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_ORIG_OUT << "\"" << html_TAB << "download=\"edata.orig.out\">Extended crystallographic data for original structure</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_RELAX_OUT << "\"" << html_TAB << "download=\"edata.relax.out\">Extended crystallographic data for relaxed structure</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_BANDS_OUT << "\"" << html_TAB << "download=\"edata.bands.out\">Extended crystallographic data for band-structure</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.agl.out"))  //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/aflow.agl.out\"" << html_TAB << "download=\"aflow.agl.out\">AGL Output</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL.out")) //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AGL.out\"" << html_TAB << "download=\"AGL.out\">AGL complete output</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_energies_temperature.out")) //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AGL_energies_temperature.out\"" << html_TAB << "download=\"AGL_energies_temperature.out\">AGL Energy versus Temperature</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_thermal_properties_temperature.out")) //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AGL_thermal_properties_temperature.out\"" << html_TAB << "download=\"AGL_thermal_properties_temperature.out\">AGL Thermal Properties versus Temperature</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"Hugoniot.out")) //CT20181212
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AGL_Hugoniot.out\"" << html_TAB << "download=\"AGL_Hugoniot.out\">AGL Hugoniot Relation</a></li>" << endl; //CT20181212
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.ael.out")) //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/aflow.ael.out\"" << html_TAB << "download=\"aflow.ael.out\">AEL Output</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Elastic_constants.out")) //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AEL_Elastic_constants.out\"" << html_TAB << "download=\"AEL_Elastic_constants.out\">AEL Elastic constants</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Compliance_tensor.out")) //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AEL_Compliance_tensor.out\"" << html_TAB << "download=\"AEL_Compliance_tensor.out\">AEL Compliance constants</a></li>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      string abader_out=aentry.prototype+"_abader.out"; //JPO20180809
//[SC20200327 - OBSOLETE]      if(aurostd::substring2bool(aentry.vfiles_WEB,abader_out)) { //JPO20180809
//[SC20200327 - OBSOLETE]        oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << abader_out << "\"" << html_TAB << "download=\"abader.out\">Bader Output</a></li>" << endl; } //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "</ul></div></div>" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]      oss << "<!-- Downloadable Files: END -->" << endl; //JPO20180809
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      // ***************************************************************************
//[SC20200327 - OBSOLETE]      // OTHER
//[SC20200327 - OBSOLETE]      if(0) if((XHOST.hostname=="nietzsche.mems.duke.edu" || XHOST.hostname=="materials.duke.edu" )&& !directory.empty()) {
//[SC20200327 - OBSOLETE]        oss << line_rule << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- DEBUG: BEGIN -->" << endl;
//[SC20200327 - OBSOLETE]        oss << "<b>DEBUG: only in " << XHOST.hostname << "</b><br>" << endl;
//[SC20200327 - OBSOLETE]        oss << aentry.entry << "<br><br>" << endl;
//[SC20200327 - OBSOLETE]        oss << line_rule << endl;
//[SC20200327 - OBSOLETE]        oss << "<!-- DEBUG: END -->" << endl;
//[SC20200327 - OBSOLETE]      }
//[SC20200327 - OBSOLETE]      // DONE
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]      oss << "<! HS WORK AFTER HERE> " << endl;
//[SC20200327 - OBSOLETE]    }
//[SC20200327 - OBSOLETE]    return voptions.size();
//[SC20200327 - OBSOLETE]  }
//[SC20200327 - OBSOLETE]}
//[SC20200327 - OBSOLETE]
//[SC20200327 - OBSOLETE]//#include "aflowlib_web_interface_test2.cpp"
//[SC20200327 - OBSOLETE]//#include "aflowlib_web_interface_test3.cpp"

// ***************************************************************************
namespace aflowlib {
  uint WEB_Aflowlib_Entry(string options,ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "aflowlib::WEB_Aflowlib_Entry():";
    if(LDEBUG) cout << XPID << "aflowlib::WEB_Aflowlib_Entry: begin<br>" << endl;

    stringstream num_prec;
    vector<string> voptions;
    aurostd::string2tokens(options,voptions,",");
    if(voptions.size()==0) {
      init::ErrorOption(options,"aflowlib::WEB_Aflowlib_Entry","aflow --aflowlib=entry");  //CO20200624 - soft patch for FR+web
      return 0; //CO20200624 - 0 is error here
    } 

    // move on
    for(uint ioption=0;ioption<voptions.size();ioption++) {
      //  oss << voptions.at(ioption) << endl;
      string option=voptions.at(ioption); //aurostd::args2attachedstring(argv,"--aflowlib=",(string) "nan");
      if(option.at(option.size()-1)=='/'|| option.at(option.size()-1)=='.') option.erase(option.end()-1,option.end()-0); //  some demoronization
      if(option.at(0)=='/'|| option.at(0)=='.') option.erase(option.begin(),option.begin()+1); //  some demoronization
      string directory="";
      string directory_RAW="",directory_LIB="",directory_WEB="";
      string directory_AUID_LIB="",directory_AUID_RAW="",directory_AUID_WEB="";
      string url_WEB;
      string label="";
      // string line_gif="<br><img border=0 width=60% height=2 src=http://materials.duke.edu/auro/images/line.gif><br><br>";
      // string art058_link=" https://doi.org/10.1016/j.commatsci.2010.05.010";
      // string art064_link=" https://doi.org/10.1021/co200012w";
      // string icsd_link=" https://www.fiz-karlsruhe.com/icsd.html";
      // string aflow_ael_readme=" http://materials.duke.edu/AFLOW/README_AFLOW_AEL.TXT"; //CO20180817
      // string art096_link=" https://doi.org/10.1103/PhysRevB.90.174107";
      // string art100_link=" https://www.nature.com/articles/sdata20159";
      // string aflow_agl_readme=" http://materials.duke.edu/AFLOW/README_AFLOW_AGL.TXT"; //CO20180817
      // string art115_link=" https://doi.org/10.1103/PhysRevMaterials.1.015401"; //CO20180817
      // string aflow_sym_readme=" http://materials.duke.edu/AFLOW/README_AFLOW_SYM.TXT"; //CO20180817
      // string art135_link=" https://doi.org/10.1107/S2053273318003066"; //CO20180817

      //DX20180817
      // string bravais_lattice_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_orig";
      // string bravais_lattice_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_relax";
      // string lattice_system_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_orig";
      // string lattice_variation_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_orig";
      // string lattice_system_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_relax";
      // string lattice_variation_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_relax";
      // string Pearson_symbol_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_orig";
      // string Pearson_symbol_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_relax";
      // string sg_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg";
      // string sg2_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg2";
      // string spacegroup_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_orig";
      // string spacegroup_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_relax";

      aflowlib::_aflowlib_entry aentry;

      xoption vflags;
      vflags.flag("FLAG::PREAMBLE",TRUE);
      vflags.flag("FLAG::CALCULATION",TRUE);
      vflags.flag("FLAG::JMOL",TRUE);
      vflags.flag("FLAG::EDATA_ORIG",FALSE);
      vflags.flag("FLAG::EDATA_RELAX",TRUE);
      vflags.flag("FLAG::THERMODYNAMICS",TRUE);
      vflags.flag("FLAG::MAGNETIC",TRUE);
      vflags.flag("FLAG::ELECTRONIC",FALSE);     // will setup later
      vflags.flag("FLAG::SCINTILLATION",TRUE);   // will setup later
      vflags.flag("FLAG::AGL",FALSE);            // will setup later
      vflags.flag("FLAG::AEL",FALSE);            // will setup later
      vflags.flag("FLAG::BADER",FALSE);          // will setup later

      // check if ICSD inside (anyway)
      string lattices[]={"BCC","BCT","CUB","FCC","HEX","MCL","MCLC","ORC","ORCC","ORCF","ORCI","RHL","TET","TRI"};
      vector<string> vline,tokens;

      string html_TAB=" target=\"_blank\"";

      if(LDEBUG) cout << "WEB_Aflowlib_Entry: [1]<br>" << endl;
      if(LDEBUG) cout << "WEB_Aflowlib_Entry: [4]<br>" << endl;
      if(LDEBUG) cout << "WEB_Aflowlib_Entry: option=" << option << endl;

      // START SEARCH

      vflags.flag("FLAG::FOUND",FALSE);
      string catalog="",auid="",strtmp,errormsg="";

      // **********************************************************************************************************
      // TRY AUID
      // **********************************************************************************************************
      if(!vflags.flag("FLAG::FOUND") && aurostd::substring2bool(aurostd::tolower(option),"aflow:")) { // CHECK AUID
        if(LDEBUG) cout << "WEB_Aflowlib_Entry: option=" << option << endl;
        string auid=aurostd::tolower(option);
        if(auid.size()!=22) {
          errormsg="aflowlib::WEB_Aflowlib_Entry:_error_on_size_of_auid="+auid;
        } else {
          // NEW
          directory_AUID_LIB=init::AFLOW_Projects_Directories("AUID")+"/"+auid.substr(0,8); for(uint i=8;i<=20;i+=2) directory_AUID_LIB+="/"+auid.substr(i,2);  // splitting aflow:ab/cd..
          directory_AUID_WEB=directory_AUID_LIB+"/WEB";
          directory_AUID_RAW=directory_AUID_LIB+"/RAW";
          directory_AUID_LIB=directory_AUID_LIB+"/LIB";
          if(LDEBUG) cout << "WEB_Aflowlib_Entry: directory_AUID_LIB=" << directory_AUID_LIB << endl;
          if(LDEBUG) cout << "WEB_Aflowlib_Entry: directory_AUID_RAW=" << directory_AUID_RAW << endl;
          if(LDEBUG) cout << "WEB_Aflowlib_Entry: directory_AUID_WEB=" << directory_AUID_WEB << endl;
          directory="";
          if(!aurostd::FileExist(directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            errormsg="aflowlib::WEB_Aflowlib_Entry:_entry_does_not_exist="+directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT;
          } else {
            _aflowlib_entry entry_tmp(string(directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT));
            directory=entry_tmp.aurl;
            auid=entry_tmp.aurl;
            aurostd::StringSubst(directory,"aflowlib.duke.edu:","");
            aurostd::StringSubst(directory,"materials.duke.edu:","");
            aurostd::StringSubst(directory,"AFLOWDATA/ICSD_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/ICSD_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB0_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB0_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB1_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB1_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB2_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB2_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB3_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB3_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB4_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB4_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB5_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB5_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB6_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB6_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB7_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB7_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB8_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB8_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB9_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB9_WEB/","");
            vflags.flag("FLAG::AUID",TRUE);
            vflags.flag("FLAG::FOUND",TRUE);
            catalog=entry_tmp.catalog;
            label=directory;
            for(uint ilat=0;ilat<14;ilat++)
              aurostd::StringSubst(label,lattices[ilat]+"/","");
            aurostd::StringSubst(label,"/",".");
          }
        }
      }

      // **********************************************************************************************************
      // TRY DIRECTORY
      // **********************************************************************************************************
      if(!vflags.flag("FLAG::FOUND")) { // tests with proto name
        string dir2test;
        dir2test=option;

        vector<string> vdir2test;
        for(uint i0=0;i0<dir2test.size();i0++) { // allow all identical indices...
          for(uint i1=i0;i1<dir2test.size();i1++) { // allow all identical indices
            for(uint i2=i1;i2<dir2test.size();i2++) { // allow all identical indices
              //	    for(uint i3=i2;i3<dir2test.size();i3++) { // allow all identical indices
              dir2test=option;
              if(dir2test.at(i0)=='.') dir2test.at(i0)='/';
              if(dir2test.at(i1)=='.') dir2test.at(i1)='/';
              if(dir2test.at(i2)=='.') dir2test.at(i2)='/';
              //      if(dir2test.at(i3)=='.') dir2test.at(i3)='/';
              bool found=false;
              for(uint i=0;i<vdir2test.size()&&!found;i++) if(dir2test==vdir2test.at(i)) found=true;
              if(!found) {
                vdir2test.push_back(dir2test);
              }
            }
            //	  }
          }
        }   
        if(LDEBUG) cout << "WEB_Aflowlib_Entry: testing dir2test=" << dir2test << endl;
        for(uint i=0;i<vdir2test.size()&&!vflags.flag("FLAG::FOUND");i++) { // allow i=j=k so that 1 OR 2 OR 3 dots are also tested
          dir2test=vdir2test.at(i);
          //		cout << "testing(" << i << ") = " << dir2test << endl;
          for(uint ilat=0;ilat<14&&!vflags.flag("FLAG::FOUND");ilat++) {
            if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("ICSD")+"/RAW/"+lattices[ilat]+"/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
              catalog="ICSD";directory=lattices[ilat]+"/"+dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);
            }
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB0")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB0";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB1")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB1";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB2";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB3")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB3";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB4")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB4";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB5")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB5";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB6")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB6";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB7")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB7";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB8")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB8";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB9")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB9";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
        }
      }

      // **********************************************************************************************************
      // TRY ICSD LINK
      // **********************************************************************************************************
      if(!vflags.flag("FLAG::FOUND")) { // icsd link
        string directory_ICSD2LINK=init::AFLOW_Projects_Directories("AUID")+"/icsd:/"+option;
        aurostd::StringSubst(directory_ICSD2LINK,"ICSD:","icsd:");
        aurostd::StringSubst(directory_ICSD2LINK,"icsd:icsd:","icsd:");    
        //	cerr << directory_ICSD2LINK << endl;
        if(aurostd::FileExist(directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
          _aflowlib_entry entry_tmp(string(directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT));
          directory=entry_tmp.aurl;
          auid=entry_tmp.aurl;
          aurostd::StringSubst(directory,"aflowlib.duke.edu:","");
          aurostd::StringSubst(directory,"materials.duke.edu:","");
          aurostd::StringSubst(directory,"AFLOWDATA/ICSD_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/ICSD_WEB/","");
          vflags.flag("FLAG::ICSD",TRUE);
          vflags.flag("FLAG::FOUND",TRUE);
          catalog=entry_tmp.catalog;
          label=directory;
          //	  cerr << directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << endl;
        }
      }

      // **********************************************************************************************************
      // CONSIDERED FOUND
      // **********************************************************************************************************    
      if(vflags.flag("FLAG::FOUND")) {
        if(catalog=="ICSD") {
          vflags.flag("FLAG::ICSD",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("ICSD")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("ICSD")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("ICSD")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/ICSD_WEB/"+directory;
        }

        if(catalog=="LIB0") {
          vflags.flag("FLAG::LIB0",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB0")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB0")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB0")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB0_RAW/"+directory;
        }
        if(vflags.flag("FLAG::FOUND") && catalog=="LIB1") {
          vflags.flag("FLAG::LIB1",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB1")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB1")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB1")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB1_RAW/"+directory;
        }
        if(catalog=="LIB2") {
          vflags.flag("FLAG::LIB2",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB2")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+directory; // June 2016
          directory_RAW=init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB2_RAW/"+directory; // May 2014
        }
        if(catalog=="LIB3") {
          vflags.flag("FLAG::LIB3",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB3")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB3")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB3")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB3_WEB/"+directory;
        }
        if(catalog=="LIB4") {
          vflags.flag("FLAG::LIB4",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB4")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB4")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB4")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB4_WEB/"+directory;
        }
        if(catalog=="LIB5") {
          vflags.flag("FLAG::LIB5",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB5")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB5")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB5")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB5_WEB/"+directory;
        }
        if(catalog=="LIB6") {
          vflags.flag("FLAG::LIB6",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB6")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB6")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB6")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB6_WEB/"+directory;
        }
        if(catalog=="LIB7") {
          vflags.flag("FLAG::LIB7",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB7")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB7")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB7")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB7_WEB/"+directory;
        }
        if(catalog=="LIB8") {
          vflags.flag("FLAG::LIB8",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB8")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB8")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB8")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB8_WEB/"+directory;
        }
        if(catalog=="LIB9") {
          vflags.flag("FLAG::LIB9",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB9")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB9")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB9")+"/RAW/"+directory;	
          url_WEB="/AFLOWDATA/LIB9_WEB/"+directory;
        }
      } else {
        errormsg="aflowlib::WEB_Aflowlib_Entry:_entry_does_not_exist="+option;
      }

      // got it  ?
      if(!aurostd::FileExist(directory_RAW+"/"+_AFLOWIN_)) directory_RAW="";

      if(!directory.empty()) { // play with aentry.entry
        aurostd::StringSubst(label,"/",".");
        aentry.file2aflowlib(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,oss);  //   oss << aentry.entry << endl;
        directory_AUID_LIB=init::AFLOW_Projects_Directories("AUID")+"/"+aflowlib::auid2directory(aentry.auid);
        directory_AUID_WEB=directory_AUID_LIB+"/WEB";
        directory_AUID_RAW=directory_AUID_LIB+"/RAW";
        directory_AUID_LIB=directory_AUID_LIB+"/LIB";    
        aurostd::string2tokens(aentry.sg2,tokens,"#");if(tokens.size()>0) aentry.sg2=tokens.at(tokens.size()-1);
        if(aentry.vfiles_WEB.size()==0) aentry.vfiles_WEB=aentry.vfiles;
      }

      if(vflags.flag("FLAG::FOUND")) {
        // check AGL/AEL
        vflags.flag("FLAG::ELECTRONIC",aurostd::substring2bool(aentry.vloop,"bands"));
        vflags.flag("FLAG::SCINTILLATION",aurostd::substring2bool(aentry.vloop,"bands"));
        vflags.flag("FLAG::AGL",aurostd::substring2bool(aentry.vloop,"agl"));
        vflags.flag("FLAG::AEL",aurostd::substring2bool(aentry.vloop,"ael"));
        vflags.flag("FLAG::BADER",aurostd::substring2bool(aentry.vloop,"bader"));
      }
      stringstream aflowlib_json;
      if(vflags.flag("FLAG::FOUND")) {
        strtmp=aurostd::efile2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON);
        aurostd::StringSubst(strtmp,"}\n",""); // remove trailing bracket add it at the end
        aflowlib_json << strtmp;
      } else {
        aflowlib_json << "{DUMMY"; // will remove at the end
      }
      stringstream aflowlib_out;
      if(vflags.flag("FLAG::FOUND")) {
        strtmp=aurostd::efile2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
        aurostd::StringSubst(strtmp,"\n",""); // remove trailing bracket add it at the end
        aflowlib_out << strtmp;
      } else {
        aflowlib_out << "DUMMY"; // will remove at the end
      }

      // adding pieces 
      //   aurostd::StringSubst(aflowlib_json,"}\n",""); // remove trailing bracket add it at the end

      // XHOST.hostname
      aflowlib_json << "," << "\"XHOST.hostname\":" << "\"" << XHOST.hostname << "\"";
      aflowlib_out << " | " << "XHOST.hostname=" << XHOST.hostname;
      // option
      aflowlib_json << "," << "\"XHOST.option\":" << "\"" << option << "\"";
      aflowlib_out << " | " << "XHOST.option=" << option;
      // label
      aflowlib_json << "," << "\"XHOST.label\":" << "\"" << label << "\"";
      aflowlib_out << " | " << "XHOST.label=" << label;
      // directory
      aflowlib_json << "," << "\"XHOST.directory\":" << "\"" << directory << "\"";
      aflowlib_out << " | " << "XHOST.directory=" << directory;

      if(vflags.flag("FLAG::FOUND")) {
        // directory_LIB
        aflowlib_json << "," << "\"XHOST.directory_LIB\":" << "\"" << directory_LIB << "\"";
        aflowlib_out << " | " << "XHOST.directory_LIB=" << directory_LIB;
        // directory_RAW
        aflowlib_json << "," << "\"XHOST.directory_RAW\":" << "\"" << directory_RAW << "\"";
        aflowlib_out << " | " << "XHOST.directory_RAW=" << directory_RAW;
        // directory_WEB
        aflowlib_json << "," << "\"XHOST.directory_WEB\":" << "\"" << directory_WEB << "\"";
        aflowlib_out << " | " << "XHOST.directory_WEB=" << directory_WEB;
        // directory_AUID_LIB
        aflowlib_json << "," << "\"XHOST.directory_AUID_LIB\":" << "\"" << directory_AUID_LIB << "\"";
        aflowlib_out << " | " << "XHOST.directory_AUID_LIB=" << directory_AUID_LIB;
        // directory_AUID_RAW
        aflowlib_json << "," << "\"XHOST.directory_AUID_RAW\":" << "\"" << directory_AUID_RAW << "\"";
        aflowlib_out << " | " << "XHOST.directory_AUID_RAW=" << directory_AUID_RAW;
        // directory_AUID_WEB
        aflowlib_json << "," << "\"XHOST.directory_AUID_WEB\":" << "\"" << directory_AUID_WEB << "\"";
        aflowlib_out << " | " << "XHOST.directory_AUID_WEB=" << directory_AUID_WEB;
      }

      if(vflags.flag("FLAG::FOUND")) {
        // XHOST.FLAG::ICSD
        aflowlib_json << "," << "\"XHOST.FLAG::ICSD\":" << (vflags.flag("FLAG::ICSD")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::ICSD=" << (vflags.flag("FLAG::ICSD")?"1":"0");
        // XHOST.FLAG::LIB0
        aflowlib_json << "," << "\"XHOST.FLAG::LIB0\":" << (vflags.flag("FLAG::LIB0")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB0=" << (vflags.flag("FLAG::LIB0")?"1":"0");
        // XHOST.FLAG::LIB1
        aflowlib_json << "," << "\"XHOST.FLAG::LIB1\":" << (vflags.flag("FLAG::LIB1")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB1=" << (vflags.flag("FLAG::LIB1")?"1":"0");
        // XHOST.FLAG::LIB2
        aflowlib_json << "," << "\"XHOST.FLAG::LIB2\":" << (vflags.flag("FLAG::LIB2")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB2=" << (vflags.flag("FLAG::LIB2")?"1":"0");
        // XHOST.FLAG::LIB3
        aflowlib_json << "," << "\"XHOST.FLAG::LIB3\":" << (vflags.flag("FLAG::LIB3")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB3=" << (vflags.flag("FLAG::LIB3")?"1":"0");
        // XHOST.FLAG::LIB4
        aflowlib_json << "," << "\"XHOST.FLAG::LIB4\":" << (vflags.flag("FLAG::LIB4")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB4=" << (vflags.flag("FLAG::LIB4")?"1":"0");
        // XHOST.FLAG::LIB5
        aflowlib_json << "," << "\"XHOST.FLAG::LIB5\":" << (vflags.flag("FLAG::LIB5")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB5=" << (vflags.flag("FLAG::LIB5")?"1":"0");
        // XHOST.FLAG::LIB6
        aflowlib_json << "," << "\"XHOST.FLAG::LIB6\":" << (vflags.flag("FLAG::LIB6")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB6=" << (vflags.flag("FLAG::LIB6")?"1":"0");
        // XHOST.FLAG::LIB7
        aflowlib_json << "," << "\"XHOST.FLAG::LIB7\":" << (vflags.flag("FLAG::LIB7")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB7=" << (vflags.flag("FLAG::LIB7")?"1":"0");
        // XHOST.FLAG::LIB8
        aflowlib_json << "," << "\"XHOST.FLAG::LIB8\":" << (vflags.flag("FLAG::LIB8")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB8=" << (vflags.flag("FLAG::LIB8")?"1":"0");
        // XHOST.FLAG::LIB9
        aflowlib_json << "," << "\"XHOST.FLAG::LIB9\":" << (vflags.flag("FLAG::LIB9")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB9=" << (vflags.flag("FLAG::LIB9")?"1":"0");
        // XHOST.FLAG::AUID
        aflowlib_json << "," << "\"XHOST.FLAG::AUID\":" << (vflags.flag("FLAG::AUID")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::AUID=" << (vflags.flag("FLAG::AUID")?"1":"0");
      }
      // XHOST.FLAG::FOUND
      aflowlib_json << "," << "\"XHOST.FLAG::FOUND\":" << (vflags.flag("FLAG::FOUND")?"true":"false");
      aflowlib_out << " | " << "XHOST.FLAG::FOUND=" << (vflags.flag("FLAG::FOUND")?"1":"0");
      //   if(!vflags.flag("FLAG::FOUND"))
      {
        // errormsg
        aflowlib_json << "," << "\"XHOST.errormsg\":" << "\"" << errormsg << "\"";
        aflowlib_out << " | " << "XHOST.errormsg=" << errormsg;
      }

      if(vflags.flag("FLAG::FOUND")) {
        // XHOST.FLAG::PREAMBLE
        aflowlib_json << "," << "\"XHOST.FLAG::PREAMBLE\":" << (vflags.flag("FLAG::PREAMBLE")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::PREAMBLE=" << (vflags.flag("FLAG::PREAMBLE")?"1":"0");
        // XHOST.FLAG::CALCULATION
        aflowlib_json << "," << "\"XHOST.FLAG::CALCULATION\":" << (vflags.flag("FLAG::CALCULATION")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::CALCULATION=" << (vflags.flag("FLAG::CALCULATION")?"1":"0");
        // XHOST.FLAG::JMOL
        aflowlib_json << "," << "\"XHOST.FLAG::JMOL\":" << (vflags.flag("FLAG::JMOL")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::JMOL=" << (vflags.flag("FLAG::JMOL")?"1":"0");
        // XHOST.FLAG::EDATA_ORIG
        aflowlib_json << "," << "\"XHOST.FLAG::EDATA_ORIG\":" << (vflags.flag("FLAG::EDATA_ORIG")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::EDATA_ORIG=" << (vflags.flag("FLAG::EDATA_ORIG")?"1":"0");
        // XHOST.FLAG::EDATA_RELAX
        aflowlib_json << "," << "\"XHOST.FLAG::EDATA_RELAX\":" << (vflags.flag("FLAG::EDATA_RELAX")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::EDATA_RELAX=" << (vflags.flag("FLAG::EDATA_RELAX")?"1":"0");
        // XHOST.FLAG::THERMODYNAMICS
        aflowlib_json << "," << "\"XHOST.FLAG::THERMODYNAMICS\":" << (vflags.flag("FLAG::THERMODYNAMICS")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::THERMODYNAMICS=" << (vflags.flag("FLAG::THERMODYNAMICS")?"1":"0");
        // XHOST.FLAG::MAGNETIC
        aflowlib_json << "," << "\"XHOST.FLAG::MAGNETIC\":" << (vflags.flag("FLAG::MAGNETIC")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::MAGNETIC=" << (vflags.flag("FLAG::MAGNETIC")?"1":"0");
        // XHOST.FLAG::ELECTRONIC
        aflowlib_json << "," << "\"XHOST.FLAG::ELECTRONIC\":" << (vflags.flag("FLAG::ELECTRONIC")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::ELECTRONIC=" << (vflags.flag("FLAG::ELECTRONIC")?"1":"0");
        // XHOST.FLAG::SCINTILLATION
        aflowlib_json << "," << "\"XHOST.FLAG::SCINTILLATION\":" << (vflags.flag("FLAG::SCINTILLATION")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::SCINTILLATION=" << (vflags.flag("FLAG::SCINTILLATION")?"1":"0");
        // XHOST.FLAG::AGL
        aflowlib_json << "," << "\"XHOST.FLAG::AGL\":" << (vflags.flag("FLAG::AGL")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::AGL=" << (vflags.flag("FLAG::AGL")?"1":"0");
        // XHOST.FLAG::AEL
        aflowlib_json << "," << "\"XHOST.FLAG::AEL\":" << (vflags.flag("FLAG::AEL")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::AEL=" << (vflags.flag("FLAG::AEL")?"1":"0");
        // XHOST.FLAG::BADER
        aflowlib_json << "," << "\"XHOST.FLAG::BADER\":" << (vflags.flag("FLAG::BADER")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::BADER=" << (vflags.flag("FLAG::BADER")?"1":"0");

        //ME20191004 START
        // Grab compressed files
        if(XHOST.vflag_control.flag("PRINT_MODE::JSON") || !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
          string content;
          // fgroup for JMOL applet
          if (aurostd::EFileExist(directory_RAW + "/aflow.fgroup.bands.json")) {
            content = aurostd::efile2string(directory_RAW + "/aflow.fgroup.bands.json");
          } else if (aurostd::EFileExist(directory_RAW + "/aflow.fgroup.relax.json")) {
            content = aurostd::efile2string(directory_RAW + "/aflow.fgroup.relax.json");
          }
          aflowlib_json << ", \"fgroup\":" << (content.empty()?"null":content);

          content = "";
          if (vflags.flag("FLAG::ELECTRONIC")) {
            //ME20200616 - Made less dependent on file name conventions
            //string system_name = KBIN::ExtractSystemName(directory_LIB);
            vector<string> vfiles;
            aurostd::DirectoryLS(directory_RAW, vfiles);
            uint nfiles = vfiles.size();
            for (uint f = 0; f < nfiles; f++) {
              if (vfiles[f].find("_bandsdata.json") != string::npos) content = aurostd::efile2string(directory_RAW + "/" + vfiles[f]);
            }
          }
          aflowlib_json << ", \"bandsdata\":" << (content.empty()?"null":content);
        }
        //ME20191004 STOP
      }

      //ME20191217 START
      // additional web output
      aflowlib_json << "," << "\"aflow_version\":\"" << AFLOW_VERSION << "\"";
      aflowlib_out << "|" << "aflow_version=" << AFLOW_VERSION;
      aflowlib_json << "," << "\"aflow_build_date\":\"" << TODAY << "\""; //CO20200624 - more semantic
      aflowlib_out << "|" << "aflow_build_date=" << TODAY;  //CO20200624 - more semantic
      //ME20191217 STOP

      // XHOST.machine_type
      aflowlib_json << "," << "\"XHOST.machine_type\":" << "\"" << XHOST.machine_type << "\"";
      aflowlib_out << " | " << "XHOST.machine_type=" << XHOST.machine_type;
      // XHOST.user
      aflowlib_json << "," << "\"XHOST.user\":" << "\"" << XHOST.user << "\"";
      aflowlib_out << " | " << "XHOST.user=" << XHOST.user;
      // XHOST.group
      aflowlib_json << "," << "\"XHOST.group\":" << "\"" << XHOST.group << "\"";
      aflowlib_out << " | " << "XHOST.group=" << XHOST.group;
      // XHOST.shell
      aflowlib_json << "," << "\"XHOST.shell\":" << "\"" << XHOST.shell << "\"";
      aflowlib_out << " | " << "XHOST.shell=" << XHOST.shell;
      // XHOST.progname
      aflowlib_json << "," << "\"XHOST.progname\":" << "\"" << XHOST.progname << "\"";
      aflowlib_out << " | " << "XHOST.progname=" << XHOST.progname;
      // XHOST.generator
      aflowlib_json << "," << "\"XHOST.generator\":" << "\"" << "aflowlib::WEB_Aflowlib_Entry" << "\"";
      aflowlib_out << " | " << "XHOST.generator=" << "aflowlib::WEB_Aflowlib_Entry";

      // wrap up
      // TXT
      if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
        oss << aflowlib_out.str() << endl;
      }
      // JSON
      if(XHOST.vflag_control.flag("PRINT_MODE::JSON") || !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
        aflowlib_json << "}";
        strtmp=aflowlib_json.str();
        aurostd::StringSubst(strtmp,"{DUMMY,","{");
        //    oss << "[" << option << "]" << endl;
        oss << strtmp << endl;
      }
    }
    return voptions.size();
  }
}


#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************

