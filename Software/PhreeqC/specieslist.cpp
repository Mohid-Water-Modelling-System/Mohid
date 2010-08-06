#include "specieslist.h"

SpeciesList::SpeciesList()
{
}

SpeciesList::~SpeciesList()
{
}

bool SpeciesList::Sort(Species *s_hplus, int start)
{
	bool swap;
	SpeciesInfo *a, *b;
	int i;

	this->s_hplus = s_hplus;

	start++;
	if (start == 0) 
		return false;

	if (start >= _count)
		return true;

	do
	{
		swap = false;
		for (i = start; i < _count; i++)
		{
			a = _index[i-1]->data;
			b = _index[i]->data;

			if (Compare(a, b) == 1) //1 means that a is greater than b
			{
				Swap(i-1, i);
				swap = true;
			}
		}
	}
	while(swap);

	return true;
}

bool SpeciesList::SortByMasterOnly(Species *s_hplus, int start)
{
	bool swap;
	SpeciesInfo *a, *b;
	int i;

	this->s_hplus = s_hplus;

	start++;
	if (start <= 0)
		return false;

	if (start >= _count)
		return true;

	do
	{
		swap = false;
		for (i = start; i < _count; i++)
		{
			a = _index[i-1]->data;
			b = _index[i]->data;

			if (CompareByMasterOnly(a, b) == 1) //1 means that a is greater than b
			{
				Swap(i-1, i);
				swap = true;
			}
		}
	}
	while(swap);

	return true;
}

int SpeciesList::Compare(SpeciesInfo *a, SpeciesInfo *b)
{
	// Put H+ first
	if (a->master_s != b->master_s)
	{
		if (a->master_s == s_hplus)
			return -1;
		if (b->master_s == s_hplus)
			return 1;
	}

  String name1, name2;

	// Other element valence states
	if (a->master_s->secondary != NULL)
		name1 = a->master_s->secondary->e->name;
	else
		name1 = a->master_s->primary->e->name;

	if (b->master_s->secondary != NULL)
		name2 = b->master_s->secondary->e->name;
	else
		name2 = b->master_s->primary->e->name;

	// Compare name of primary or secondary master species; log molality
	int j = name2.Compare(name1);

	if (j != 0) 
		return j;

	if (a->s->lm > b->s->lm)
		return -1;
	else if (a->s->lm < b->s->lm)
		return 1;
	else
		return 0;
}

int SpeciesList::CompareByMasterOnly(SpeciesInfo *a, SpeciesInfo *b)
{
	// Put H+ first
	if (a->master_s != b->master_s)
	{
		if (a->master_s == s_hplus)
			return -1;
		if (b->master_s == s_hplus)
			return 1;
	}

  String name1, name2;

	// Other element valence states
	if (a->master_s->secondary != NULL)
		name1 = a->master_s->secondary->e->name;
	else
		name1 = a->master_s->primary->e->name;

	if (b->master_s->secondary != NULL)
		name2 = b->master_s->secondary->e->name;
	else
		name2 = b->master_s->primary->e->name;

	// Compare name of primary or secondary master species; log molality
	return name2.Compare(name1);
}