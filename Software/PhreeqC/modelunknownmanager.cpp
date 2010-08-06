#include "modelunknownmanager.h"
#include <string.h>

ModelUnknownManager::ModelUnknownManager()
{
	md = NULL;
	//unknown_list_ptr = NULL;
	unknown_list = NULL;

	count_unknowns = 0;

	unknown_list = new List<Unknown>;
}

ModelUnknownManager::~ModelUnknownManager()
{
	delete unknown_list;
}

void ModelUnknownManager::SetData(ModelData *md_ptr)
{
	md = md_ptr;
}

void ModelUnknownManager::SetupUnknowns()
{
	//md MUST me defined using SetData

	int i, j;
	Use *use;

  md->max_unknowns = 0;

	use = &md->use;

	// Count mass balance in solution
	md->max_unknowns += use->sol_p->totals->Count();

	// Add 5 for ionic strength, activity of water, charge balance, total H, total O
  md->max_unknowns += 5;

	// Count pure phases
	if (use->ppa_p != NULL)
		md->max_unknowns += use->ppa_p->pure_phases->Count();

	// Count exchange
	if (md->use.exc_p != NULL)
	{
		ExchComp *excc_p;
		ElementOfSpecies *eos_p;

		for (i = 0; i < use->exc_p->comps->Count(); i++)
		{
			excc_p = (*use->exc_p->comps)[i];

			for (j = 0; j < excc_p->totals->Count(); j++)
			{
				eos_p = (*excc_p->totals)[j];

				if (eos_p->e->master == NULL) //<= maybe this verification isn't necessary
				{
					sprintf(message, "Master species missing for element %s", eos_p->e->name.CharPtr());
					throw ExceptionHandler(message);
				}

				if (eos_p->e->master->type == EX)
					md->max_unknowns++;
			}
		}
	}

	// Count surfaces
  if (use->sur_p != NULL)
  {
    if (use->sur_p->type != CD_MUSIC)
			md->max_unknowns += use->sur_p->comps->Count() + use->sur_p->charge->Count();
    else
			md->max_unknowns +=	use->sur_p->comps->Count() + (4 * use->sur_p->charge->Count());
  }

	// Count gas components
	if (use->gas_p != NULL)
		md->max_unknowns++;

	// Count solid solutions
	if (use->ssa_p != NULL)
	{
		SS *ss_p;

		for (i = 0; i < use->ssa_p->ss_list->Count(); i++)
		{
			ss_p = (*use->ssa_p->ss_list)[i];
			md->max_unknowns += ss_p->comps_list->Count();
		}
	}

	//One for luck
  md->max_unknowns++;

	// "Allocate" space for pointer array and structures
	unknown_list->Clear();
	unknown_list->SetNewCapacity(md->max_unknowns);

	Unknown *u;
	for (i = 0; i < md->max_unknowns; i++)
	{
		u = unknown_list->AddNew(); //ToDo: Check Unknown Constructor with "unknown_alloc" original function to see if it's all ok
		u->number = i;
	}

	//unknown_list_ptr = unknown_list->Pointer();
}