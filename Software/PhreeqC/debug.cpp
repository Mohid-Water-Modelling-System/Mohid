#include "debug.h"
#include <math.h>

Debug::Debug()
{
	debug_status = false;

	addsolution_f = NULL;
	addexchange_f = NULL;
	addppassemblage_f = NULL;
}

Debug::~Debug()
{
}

void Debug::OpenSpeciesReactionFile(String filename)
{
	srxn_f = fopen(filename.CharPtr(), "wt");
}

void Debug::CloseSpeciesReactionFile()
{
	fclose(srxn_f);
}

void Debug::PrintSpeciesReaction(Species *s)
{
	ReactionToken *rt;

	fprintf(srxn_f, "Species: %s\n", s->name.CharPtr());
	for (int i = 0; i < s->rxn->token_list->Count(); i++)
	{
		rt = s->rxn->token_list->Element(i);

		fprintf(srxn_f, "token: %s, coef: %f\n", rt->name.CharPtr(), rt->coef);
	}
}

void Debug::PrintMasterList(String filename, GlobalData *gd)
{
	FILE *f = fopen (filename.CharPtr(), "wt");

	for (int i = 0; i < gd->master_list.Count(); i++)
	{
		fprintf(f, "%d:", i);
		gd->master_list[i]->PrintToFile(f, 1);
		fprintf(f, "=====================================\n");
	}

	fclose(f);
}

void Debug::PrintSpeciesList(String filename, GlobalData *gd)
{
	FILE *f = fopen (filename.CharPtr(), "wt");

	for (int i = 0; i < gd->species_list.Count(); i++)
	{
		fprintf(f, "%d: ", i);
		gd->species_list[i]->PrintToFile(f, 1);
		fprintf(f, "=====================================\n");
	}

	fclose(f);
}

void Debug::PrintUnknownInfoList(String filename, List<UnknownInfo> *ui)
{
	FILE *f = fopen (filename.CharPtr(), "wt");

	for (int i = 0; i < ui->Count(); i++)
	{
		fprintf(f, "%d:\n", i);
		(*ui)[i]->PrintToFile(f, 1);
		fprintf(f, "=====================================\n");
	}

	fclose(f);
}

void Debug::PrintArray(String filename, LDBLE *array_, int count)
{
	FILE *f = fopen (filename.CharPtr(), "wt");

	for (int i = 0; i < count; i++)
		fprintf(f, "%03d: %f\n", i, array_[i]);

	fclose(f);
}

void Debug::PrintUnknownList(String filename, GlobalData *gd)
{
	FILE *f = fopen (filename.CharPtr(), "wt");

	/*
	for (int i = 0; i < gd->x.Count(); i++)
	{
		gd->x[i]->PrintToFile(f, 1);
		fprintf(f, "=====================================\n");
	}
	*/

	fclose(f);
}

void Debug::OpenUnknownFile(String filename)
{
	u_f = fopen(filename.CharPtr(), "wt");
}

void Debug::CloseUnknownFile()
{
	fclose(u_f);
}

void Debug::PrintUnknownToFile(List<Unknown> *x, String comment)
{
	fprintf(u_f, "\nINICIO %s================================\n\n", comment.CharPtr());
	for (int i = 0; i < x->Count(); i++)
	{
		fprintf(u_f, "%d: \n", i);
		(*x)[i]->PrintToFile(u_f, 1);
		fprintf(u_f, "-------------------------------------\n");
	}
	fprintf(u_f, "\nFIM %s================================\n", comment.CharPtr());
}

void Debug::PrintPhaseList(String filename, List<Phase> *p)
{
	FILE *f = fopen (filename.CharPtr(), "wt");

	for (int i = 0; i < p->Count(); i++)
	{
		fprintf(f, "%d:", i);
		(*p)[i]->PrintToFile(f, 1);
		fprintf(f, "=====================================\n");
	}

	fclose(f);
}

void Debug::OpenArrayFile(String filename)
{
	arr_f = fopen(filename.CharPtr(), "wt");
	arr_count = 0;
}

void Debug::CloseArrayFile()
{
	fclose(arr_f);
}

void Debug::PrintArrayToFile(LDBLE *array_, int count, String comment)
{
	fprintf(arr_f, "\nINICIO %03d: %s================================\n\n", arr_count, comment.CharPtr());

	arr_count++;
	for (int i = 0; i < count; i++)
		fprintf(arr_f, "%03d: %f\n", i, array_[i]);

	fprintf(arr_f, "\nFIM %s================================\n", comment.CharPtr());
}

void Debug::PrintSpeciesInfoListToFile(String filename, SpeciesList *sl)
{
	FILE *f = fopen (filename.CharPtr(), "wt");

	for (int i = 0; i < sl->Count(); i++)
	{
		fprintf(f, "%d:", i);
		(*sl)[i]->s->PrintToFile(f, 1);
		fprintf(f, "=====================================\n");
	}

	fclose(f);
}

void Debug::OpenSpeciesInfoFile(String filename)
{
	si_f = fopen(filename.CharPtr(), "wt");
	si_count = 0;
}

void Debug::PrintSpeciesInfoToFile(String comment, SpeciesList *si_list, ModelData *md, GlobalData *gd)
{
	fprintf(si_f, "%s %d: %20.10e / mu_x: %20.20e\n", comment.CharPtr(), si_count, md->mb_unknown->f, md->mu_x);
	si_count++;
	/*
	fprintf(si_f, "Number: %d---------------\n%s\n---------------\n", si_count++, comment.CharPtr());

	SpeciesInfo *s;
	Master *master_ptr;
	String name1;
	String name;
	LDBLE min, lm;
	for (int i = 0; i < si_list->Count(); i++)
	{
		s = (*si_list)[i];

    if (s->s->type == EX)
      continue;
    if (s->s->type == SURF)
      continue;
    if (s->master_s->secondary != NULL)
    {
      master_ptr = s->master_s->secondary;
      name1 = s->master_s->secondary->e->name;
    }
    else
    {
      master_ptr = s->master_s->primary;
      name1 = s->master_s->primary->e->name;
    }

    if (name1 != name)
    {
      name = name1;
			fprintf(si_f, "%-14s%12.9e\n", name.CharPtr(), (double) (master_ptr->total / md->mass_water_aq_x));
      min = md->censor * master_ptr->total / md->mass_water_aq_x;
      if (min > 0)
      {
				min = log10 (min);
      }
      else
      {
				min = -1000.;
      }
    }

    if (s->s->lm > min)
    {
      if (s->s == gd->s_h2o)
      {
				lm = log10 (gd->s_h2o->moles / md->mass_water_aq_x);
      }
      else
      {
				lm = s->s->lm;
      }

			fprintf(si_f, "   %-20s%12.9e  %12.9e\n", s->s->name.CharPtr(),
						  (double) ((s->s->moles) / md->mass_water_aq_x), 
							(double) lm);
    }
  }
	fprintf(si_f, "\n");	
	fprintf(si_f, "=====================================\n");
	*/
}

void Debug::CloseSpeciesInfoFile()
{
	fclose(si_f);
}

void Debug::CommentIntoSpeciesInfoFile(String comment)
{
	fprintf(si_f, "%s\n", comment.CharPtr());
}

LDBLE Debug::Under(LDBLE xval, LDBLE LOG_10)
{
	//
	// Exponentiate a number, x, but censor large and small numbers
	// log values less than MIN_LM are set to 0
	// log values greater than MAX_LM are set to 10**MAX_LM
	//

  if (xval < -40.)
    return (0.0);

	if (xval > 3.)
    return (1.0e3);

	return (exp (xval * LOG_10));
}

void Debug::OpenGammasFile()
{
	gammas_f = fopen("gammas.txt", "wt");
	gammas_count = 0;
}

void Debug::CommentIntoGammasFile(String comment)
{
	fprintf(gammas_f, "%s", comment.CharPtr());
}

void Debug::PrintGammasToFile(String comment, List<Species> *list)
{
	Species *s;
	for (int i = 0; i < list->Count(); i++)
	{
		s = (*list)[i];
		fprintf(gammas_f, "%d (%d): Species: %s / lg value: %20.20e\n", i, gammas_count, s->name.CharPtr(), s->lg);
	}
		
	gammas_count++;
	fprintf(gammas_f, "%s =====================================\n", comment.CharPtr());
}

void Debug::CloseGammasFile()
{
	fclose(gammas_f);
}

void Debug::OpenDeltaFile()
{
	delta_f = fopen("delta.txt", "wt");
	delta_count = 0;
}

void Debug::CommentIntoDeltaFile(String comment)
{
	fprintf(delta_f, "%s", comment);
}

void Debug::PrintDeltaToFile(String comment, LDBLE *delta, int count)
{
	int i;

	for (i = 0; i < count; i++)
	{
		fprintf(delta_f, "%d: %20.20e\n", delta_count, delta[i]);
	}
	delta_count++;
}

void Debug::CloseDeltaFile()
{
	fclose(delta_f);
}

void Debug::OpenLAFile()
{
	la_f = fopen("la.txt", "wt");
	la_count = 0;
}

void Debug::CommentIntoLAFile(String comment)
{
	fprintf(la_f, "%s", comment);
}

void Debug::PrintLAToFile(String comment, int count, List<Unknown> *x)
{
	int i;

	fprintf(la_f, "%s\n", comment.CharPtr());
	for (i = 0; i < count; i++)
	{
		fprintf(la_f, "%d: %20.20e\n", la_count, (*x)[i]->f);
	}
	la_count++;
}

void Debug::CloseLAFile()
{
	fclose(la_f);
}

void Debug::OpenDebugFile()
{
	debug_file = fopen("debug.txt", "wt");
}

void Debug::CloseDebugFile()
{
	fclose(debug_file);
}

void Debug::OpenSetFile()
{
	set_f = fopen("set.txt", "wt");
	set_count = 0;
}

void Debug::CloseSetFile()
{
	fclose(set_f);
}

void Debug::OpenReviseGuessesFile()
{
	reviseguesses_f = fopen("reviseguesses.txt", "wt");
	reviseguesses_count = 0;
}

void Debug::CloseReviseGuessesFile()
{
	fclose(reviseguesses_f);
}

void Debug::OpenBuildModelFile()
{
	buildmodel_f = fopen("buildmodel.txt", "wt");
	buildmodel_count = 0;
}

void Debug::CloseBuildModelFile()
{
	fclose(buildmodel_f);
}

void Debug::OpenMolalitiesFile()
{
	molalities_f = fopen("molalities.txt", "wt");
	molalities_count = 0;
}

void Debug::CloseMolalitiesFile()
{
	fclose(molalities_f);
}

void Debug::OpenMBSumsFile()
{
	mbsums_f = fopen("mbsums.txt", "wt");
	mbsums_count = 0;
}

void Debug::CloseMBSumsFile()
{
	fclose(mbsums_f);
}

void Debug::OpenBuildSSAssemblageFile()
{
	buildssassemblage_f = fopen("buildssassemblage.txt", "wt");
	buildssassemblage_count = 0;
}

void Debug::CloseBuildSSAssemblageFile()
{
	fclose(buildssassemblage_f);
}

void Debug::OpenBuildJacobianSumsFile()
{
	buildjacobiansums_f = fopen("buildjacobiansums.txt", "wt");
	buildjacobiansums_count = 0;
}

void Debug::CloseBuildJacobianSumsFile()
{
	fclose(buildjacobiansums_f);
}

void Debug::OpenIneqFile()
{
	ineq_f = fopen("ineq.txt", "wt");
	ineq_count = 0;
}

void Debug::CloseIneqFile()
{
	fclose(ineq_f);
}

void Debug::OpenResetFile()
{
	reset_f = fopen("reset.txt", "wt");
	reset_count = 0;
}

void Debug::CloseResetFile()
{
	fclose(reset_f);
}

void Debug::OpenSetupSolutionFile()
{
	setup_solution_f = fopen("setup_solution.txt", "wt");
	setup_solution_count = 0;
}

void Debug::CloseSetupSolutionFile()
{
	fclose(setup_solution_f);
}

//======================================================================================================
//
//======================================================================================================
void Debug::OpenAddSolutionFile()
{
	addsolution_f = NULL;
	addsolution_count = 0;

	addsolution_f = fopen("addsolution.txt", "wt");
}

void Debug::CloseAddSolutionFile()
{
	if (addsolution_f != NULL) fclose(addsolution_f);
}

void Debug::OpenAddExchangeFile()
{
	addexchange_f = NULL;
	addexchange_count = 0;

	addexchange_f = fopen("addexchange.txt", "wt");
}

void Debug::CloseAddExchangeFile()
{
	if (addexchange_f != NULL) fclose(addexchange_f);
}

void Debug::OpenAddPPAssemblageFile()
{
	addppassemblage_f = NULL;
	addppassemblage_count = 0;

	addexchange_f = fopen("addppassemblage.txt", "wt");
}

void Debug::CloseAddPPAssemblageFile()
{
	if (addppassemblage_f != NULL) fclose(addppassemblage_f);
}