#ifndef interfaceH
#define interfaceH

void pm_start(int *instance_id, int *status);
void pm_model_option_printalways(int *instance_id, int *print_always, int *status);
void pm_solution_data(int *instance_id, LDBLE *temperature, int *units, LDBLE *density, LDBLE *mass_water, int *status);
void pm_solution_redox(int *instance_id, char *element_1, LDBLE *valence_1, char *element_2, LDBLE *valence_2, int *status, int c1, int c2);
void pm_solution_ph(int *instance_id, LDBLE *value, int *charge, int *use_phase, char *phase_name, LDBLE *si, int *status, int c);
void pm_solution_pe(int *instance_id, LDBLE *value, int *charge, int *use_phase, char *phase_name, LDBLE *si, int *status, int c);
void pm_conc_add(int *instance_id, char *name, LDBLE *value, int *status, int c);
void pm_conc_use(int *instance_id, int *use_charge, int *use_phase, int *use_as, int *use_gfw, int *use_unit, int *use_redox, int *status);
void pm_conc_as(int *instance_id, char *as, int *status, int c);
void pm_conc_gfw(int *instance_id, LDBLE *gfw, int *status);
void pm_conc_unit(int *instance_id, int *unit, int *status);
void pm_conc_phase(int *instance_id, char *phase_name, LDBLE *si, int *status);
void pm_conc_redox(int *instance_id, char *element_1, LDBLE *valence_1, char *element_2, LDBLE *valence_2, int *status, int c1, int c2);
void pm_conc_save(int *instance_id, int *conc_id, int *ma_index, LDBLE *gfw, int *status);
void pm_use_ppa(int *instance_id, int *use, int *status);
void pm_use_exa(int *instance_id, int *use, int *status);
void pm_use_sa(int *instance_id, int *use, int *status);
void pm_set_use(int *instance_id, int *use_ppa, int *use_gas, int *use_ssa, int *use_sa, int *use_exa, int *status);
void pm_ppa_pp(int *instance_id, char *name, char *add_formula, LDBLE *si, LDBLE *moles, int *force_equality, int *dissolve_only, int *pp_id, int *status,	int c1, int c2);
void pm_exa_exchanger(int *instance_id, char *formula, int *ex_type, char *has, LDBLE *amount, int *index, int *status, int c1, int c2);
void pm_sa_options(int *instance_id, int *no_edl, int *diffuse_layer, LDBLE *diffuse_layer_thickness, int *only_counter_ions, int *status);
void pm_sa_surface(int *instance_id, int *s_type, char *name, char *that_has, LDBLE *specific_area, LDBLE *sites, LDBLE *mass, int *use_area, int *use_mass, int *status, int c1, int c2);
void pm_solution_set_temperature(int *instance_id, LDBLE *new_value, int *status);
void pm_solution_set_water_mass(int *instance_id, LDBLE *new_value, int *status);
void pm_solution_set_ph(int *instance_id, LDBLE *new_value, int *status);
void pm_solution_set_pe(int *instance_id, LDBLE *new_value, int *status);
void pm_solution_set_conc_value(int *instance_id, int *index, LDBLE *new_value, int *status);
void pm_set_exchange_value (int *instance_id, int *index, LDBLE *value, int *status);
void pm_set_phase_value (int *instance_id, int *index, LDBLE *value, int *status);
void pm_get_masterspecies_value(int *instance_id, int *index, LDBLE *value, int *MassOrConcentration, int *status);
void pm_get_species_value(int *instance_id, int *index, LDBLE *value, int *MassOrConcentration, int *status);
void pm_get_exchange_value(int *instance_id, int *index, LDBLE *value, int *status);
void pm_get_phase_value(int *instance_id, int *index, LDBLE *value, int *status);
void pm_get_species_index(int *instance_id, char *name, int *index, int *status, int c);
void pm_get_masterspecies_index(int *instance_id, char *name, int *index, int *status, int c);
void pm_run_model(int *instance_id, int *status);
void pm_run_model(int *instance_id, int *print_input, int *print_result, int *status);
void pm_run_model(int *instance_id, int *print_input, int *status);
void pm_setup_model(int *instance_id, int *status);
void pm_read_database(int *instance_id, char *database_path, int *status);
void pm_kill(int *instance_id, int *status);
void pm_set_input_value(int *instance_id, int *input_id, int *group, LDBLE *new_value, int *status);
void pm_get_result_value(int *instance_id, int *result_id, int *group, LDBLE *new_value, int *status);
void pm_get_data(int *instance_id, LDBLE *mass_of_water, LDBLE *ph, LDBLE *pe, int *status);
void pm_set_required_data(int *instance_id, LDBLE *mass_of_water, LDBLE *temp, LDBLE *density, int *status);
void pm_set_ph(int *instance_id, int *charge, LDBLE *ph, int *status);
void pm_set_pe(int *instance_id, int *charge, LDBLE *pe_value, int *status);

#endif