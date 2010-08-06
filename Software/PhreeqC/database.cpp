#include "database.h"
#include <stdio.h>
#include <conio.h>

Database::Database(GlobalData *gd): 
	DBEngine(gd)
{
	file = NULL;
	file = new File(&line);
}

Database::~Database()
{
	if (file != NULL) delete file;
}

bool Database::Load(String &file_name)
{
	string f_name; 

	f_name.assign(file_name.CharPtr()); 

	file->OpenFile(f_name);
	if (!file->IsOpen())
	{
		printf("\nO arquivo de dados %s não foi encontrado...\n", file_name.CharPtr());
		//_getch();
		return false;
	}

	//d.OpenSpeciesReactionFile("srxn.txt");

	RETURN_TYPE lt;	
	do
	{
		file->GetLine();
		if (file->Eof())
			return false;

		lt = line.LineType();
	}	
	while (lt != _keyword_);

	bool stop = false;
	while (!stop)
	{		
		switch (line.KeywordId())
		{
			case 0: //end
				stop = true;
        break;
			case 1: //solution_species
				if (!LoadSolutionSpeciesBlock())
					return false;
				//PrintToFileSpeciesInfo("species_2.txt");
				break;
			case 2: //solution_master_species
				if (!LoadSolutionMasterSpeciesBlock())
					return false;
				//PrintToFileSpeciesInfo("masterspecies_2.txt");
				break;
			case 3: //phases
				if (!LoadPhasesBlock())
					return false;
				//PrintToFileSpeciesInfo("phases_2.txt");
				break;
			case 4: //exchange_master_species
				if (!LoadExchangeMasterSpeciesBlock())
					return false;
				//PrintToFileSpeciesInfo("exchangemasterspecies_2.txt");
				break;
			case 5: //exchange_species
				if (!LoadExchangeSpeciesBlock())
					return false;
				//PrintToFileSpeciesInfo("exchangespecies_2.txt");
				break;
			case 6: //surface_master_species
				if (!LoadSurfaceMasterSpeciesBlock())
					return false;
				//PrintToFileSpeciesInfo("surfacemasterspecies_2.txt");
				break;
			case 7: //surface_species
				if (!LoadSurfaceSpeciesBlock())
					return false;
				//PrintToFileSpeciesInfo("surfacespecies_2.txt");
				break;
			default:
				stop = true;
				break;
		}

		lt = line.LineType();
		if (lt == _error_)
			return false;		
	}

	/*
	//For debug purposes
	for (int j = 0; j < gd->species_list.Count(); j++)
	{
		Species *s_ptr = gd->species_list[j];

		printf("%s: eos_list->Count(): %d\n", s_ptr->name.CharPtr(), s_ptr->eos_list->Count());

	}
	getch();
	*/


	gd->species_list.Sort(0, true);
	gd->phase_list.Sort(0,	true);

	//PrintToFilePhasesInfo("c:\\phase_beftidy_2.txt");

	//bool r = TidyDatabase();

	//gd->species_list.Sort(0, false);
	//gd->phase_list.Sort(0, false);

	//PrintToFilePhasesInfo("c:\\phase_afttidy_2.txt");

	//d.CloseSpeciesReactionFile();

	//d.PrintMasterList("master_0.txt", gd);
	//d.PrintSpeciesList("species_0.txt", gd);
	return true;
}

bool Database::LoadSolutionMasterSpeciesBlock()
{
	RETURN_TYPE lt;
	Master *m;

	for (;;)
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();

		switch (lt)
		{
		case _error_:
			return false;
		case _keyword_:
			return true;
		case _empty_:
			continue;
		}

		if (!line.NextToken() || line.TokenType() != _upper_)
			return false;

		line.Token().Replace("(+", "(");

		m = gd->master_list.Store(&line.Token());

		m->e = gd->element_list.Store(&line.Token(), false);
		m->name = m->e->name;
		
		if (!line.NextToken() || (line.TokenType() != _upper_ && line.Token().Compare("e-", true) != 0))
				return false;

		String token = line.Token();
		CorrectPlusOne(token);

		m->s = gd->species_list.Search(&token);
		if (m->s == NULL)
		{
			m->s = gd->species_list.Store(&token, false);;
			GetZ(token, m->s->z);
		}

		if (!line.NextToken())
			return false;

		m->alk = (LDBLE)line.Token().AsDouble();

		if (!line.NextToken())
			return false;
		
		if (line.TokenType() == _digit_)
			m->gfw = (LDBLE)line.Token().AsDouble();
		else
			m->gfw_formula = line.Token();

		if (m->name.Find("(") == -1)
		{
			m->primary = true;

			if (m->name.Find("E") != 0)
			{
				if (line.NextToken() && line.TokenType() == _digit_)
					m->e->gfw = (LDBLE)line.Token().AsDouble();
				else
					return false;
			}
		}
		else
			m->primary = false;
	}

	return true;
}

bool Database::LoadSolutionSpeciesBlock()
{
	options_list.Clear();
	options_list.AddNew("-no_check");
	options_list.AddNew("-check");
	options_list.AddNew("-gamma");
	options_list.AddNew("-mb");
	options_list.AddNew("-mass_balance");
	options_list.AddNew("-log_k");
	options_list.AddNew("-logk");
	options_list.AddNew("-delta_h");
	options_list.AddNew("-deltah");
	options_list.AddNew("-analytical_expression");
	options_list.AddNew("-a_e");
	options_list.AddNew("-ae");
	options_list.AddNew("-mole_balance");
	options_list.AddNew("-llnl_gamma");
	options_list.AddNew("-co2_llnl_gamma");
	options_list.AddNew("-activity_water");
	options_list.AddNew("-add_logk");
	options_list.AddNew("-add_log_k");
	options_list.AddNew("-add_constant");
	options_list.AddNew("-dw");

	s = NULL;
	RETURN_TYPE lt;

	for (;;)
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();

		switch (lt)
		{
		case _error_:
			return false;
		case _keyword_:
			return true;
		case _empty_:
			continue;
		case _option_:
			if (!LoadSolutionSpeciesOption())
				return false;
			break;
		case _none_:
			if (!LoadSolutionSpeciesEquation())
				return false;
			break;
		}
	}
	return true;
}

bool Database::LoadSolutionSpeciesOption()
{
	if (s == NULL)
		return false;

	line.NextToken();

	int option_id;
	options_list.Search(&line.Token(), option_id, true, true);

	switch(option_id)
	{
	case 0: //-no_check
		s->check_equation = false;
		break;

	case 1: //-check
		s->check_equation = true;
		break;

	case 2: //-gamma
		s->gflag = 2; //Wateq D-H
		if (!GetValueFromNextToken(s->dha) || !GetValueFromNextToken(s->dhb))
			return false;
		break;

	case 3:  //-mb
	case 4:  //-mass_balance
	case 12: //-mole_balance
		if(!ReadMassBalance(s))
			return false;
		break;

	case 5: //-log_k
	case 6: //-logk
		if (!GetValueFromNextToken(s->logk[0]))
			return false;
		break;

	case 7: //-delta_h
	case 8: //-deltah
		if (!ReadDeltaH(s->logk[1], s->original_units))
			return false;
		break;

	case 9:  //-analytical_expression
	case 10: //-a_e
	case 11: //-ae
		if (!ReadAnalyticalExpression(&(s->logk[2])))
			return false;
		break;

	case 13: //-llnl_gamma
		s->gflag = 7;		// llnl D-H
		if (!GetValueFromNextToken(s->dha))
			return false;
		break;

	case 14: //-co2_llnl_gamma
		s->gflag = 8; // llnl CO2 D-H
		break;

	case 15: //-activity_water
    s->gflag = 9;		// activity_water/55.5 for HDO, D2O, H2[O18], etc
		break;

	case 16: //-add_logk
	case 17: //-add_log_k
		if(!ReadAddLogK(s->add_logk))
			return false;
		break;

	case 18: //-add_constant
		if (!ReadAddConstant(s->add_logk))
			return false;
		break;

	case 19:
		if (!GetValueFromNextToken(s->dw))
			return false;
		break;
	}

	return true;
}

bool Database::LoadSolutionSpeciesEquation()
{
	long i;
	bool r = true;

	s = NULL;
	eos_list->Clear();
	
	if ((r = ParseEq(line.text)))
	{
		rt.token_list->Element(0)->s = gd->species_list.Store(&rt.token_list->Element(0)->name);
		rt.token_list->Element(0)->s->z = rt.token_list->Element(0)->z;
		for (i = 1; i < rt.token_list->Count(); i++)
		{
			rt.token_list->Element(i)->s = gd->species_list.Store(&rt.token_list->Element(i)->name, false);
			rt.token_list->Element(i)->s->z = rt.token_list->Element(i)->z;
		}

		s = rt.token_list->Element(0)->s;	

		eos_list->CopyTo(s->eos_list);

		ElementOfSpecies *e;
		for (i = 0; i < eos_list->Count(); i++)
		{
			e = (*eos_list)[i];
			if (e->e->name == "C")
				s->carbon = e->coef;
			else if (e->e->name == "H")
				s->h = e->coef;
			else if (e->e->name == "O")
				s->o = e->coef;
		}


		rt.CopyReactionsToSpecies(s->rxn);
		rt.Reset();

		s->dha = 0.0;
		s->dhb = 0.0;

		if (Equal(s->z, (LDBLE)0.0, (LDBLE)TOLERANCE))
		{
			s->gflag = 0;
			s->dhb   = (LDBLE)0.1;
		}
		else
			s->gflag = 1;

		if (s->name == "H+")
		{
			gd->s_hplus = s;
			gd->s_hplus->type = HPLUS;
		}
		else if (s->name == "H3O+")
		{
			gd->s_h3oplus = s;
			gd->s_h3oplus->type = HPLUS;
		}
		else if (s->name == "e-")
		{
			gd->s_eminus = s;
			gd->s_eminus->type = EMINUS;
			gd->s_eminus->gflag = 3;
		}
		else if (s->name == "H2O")
		{
			gd->s_h2o = s;
			gd->s_h2o->type = H2O;
			gd->s_h2o->gflag = 3;
		}
		else if (s->name.Find("(s)") != -1)
			s->type = SOLID;
		else if (s->name == "H2")
		{
			gd->s_h2 = s;
			gd->s_h2->type = AQ;
		}
		else if (s->name == "O2")
		{
			gd->s_o2 = s;
			gd->s_o2->type = AQ;
		}
		else
			s->type = AQ;
	}

	//d.PrintSpeciesReaction(s);

	return r;
}

bool Database::LoadPhasesBlock()
{
	options_list.Clear();
	options_list.AddNew("-no_check");
	options_list.AddNew("-check");
	options_list.AddNew("-log_k");
	options_list.AddNew("-logk");
	options_list.AddNew("-delta_h");
	options_list.AddNew("-deltah");
	options_list.AddNew("-analytical_expression");
	options_list.AddNew("-a_e");
	options_list.AddNew("-ae");
	options_list.AddNew("-add_logk");
	options_list.AddNew("-add_log_k");
	options_list.AddNew("-add_constant");

	p = NULL;

	RETURN_TYPE lt;

	for (;;)
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();

		switch (lt)
		{
		case _error_:
			return false;
		case _keyword_:
			return true;
		case _empty_:
			continue;
		case _option_:
			if (!LoadPhaseOption())
				return false;
			break;
		case _none_:
			if (!LoadPhaseEquation())
				return false;
			break;
		}
	}

	return true;
}

bool Database::LoadPhaseOption()
{
	if (p == NULL)
		return false;

	line.NextToken();

	int option_id;
	options_list.Search(&line.Token(), option_id, true, true);

	switch(option_id)
	{
	case 0: //-no_check
		p->check_equation = false;
		break;

	case 1: //-check
		p->check_equation = true;
		break;

	case 2: //-log_k
	case 3: //-logk
		if (!GetValueFromNextToken(p->logk[0]))
			return false;
		break;

	case 4: //-delta_h
	case 5: //-deltah
		if (!ReadDeltaH(p->logk[1], p->original_units))
			return false;
		break;

	case 6:  //-analytical_expression
	case 7: //-a_e
	case 8: //-ae
		if (!ReadAnalyticalExpression(&(p->logk[2])))
			return false;
		break;

	case 9: //-add_logk
	case 10: //-add_log_k
		if(!ReadAddLogK(p->add_logk))
			return false;
		break;

	case 11: //-add_constant
		if (!ReadAddConstant(p->add_logk))
			return false;
		break;

	}

	return true;
}

bool Database::LoadPhaseEquation()
{
	if (!line.NextToken())
		return false;

	p = gd->phase_list.Store(&line.Token());

	RETURN_TYPE lt;
	do
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();
	}
	while (lt == _empty_);

	if (lt != _none_)
		return false;

	eos_list->Clear();
	
	bool r;
	if ((r = ParseEq(line.text, false)))
	{
		p->formula = rt.token_list->Element(0)->name;  
		p->formula.Replace("(g)", "");
		p->formula.Replace("(s)", "");
		p->formula.Replace("(G)", "");
		p->formula.Replace("(S)", "");
		
		String s_name;
		ReactionToken *rxnt;
		for (int i = 1; i < rt.token_list->Count(); i++)
    {
			rxnt = rt.token_list->Element(i);
			if (rxnt->name.Find("(s)", true) == -1 && rxnt->name.Find("(g)", true) == -1)
			{
				s_name = rxnt->name;
				s_name.Replace("(aq)", "");
				s_name.Replace("(AQ)", "");
				s_name.Replace("H2O(l)", "H2O");
				s_name.Replace("(H2O(L)", "H2O");
				rxnt->s = gd->species_list.Store(&s_name, false);
				rxnt->s->z = rxnt->z;
			}
			else
				rxnt->s = NULL;
		}

		p->rxn->Reset();
		rt.CopyReactionsToPhase(p->rxn);

		eos_list->CopyTo(p->eos_list);
    p->type = SOLID;
	}
	else
		r = false;

	return r;
}

bool Database::LoadExchangeMasterSpeciesBlock()
{
	RETURN_TYPE lt;

	for (;;)
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();

		switch (lt)
		{
		case _error_:
		case _option_:
			return false;
		case _keyword_:
			return true;
		case _empty_:
			continue;
		}

		if (!line.NextToken())
			return false;

		line.Token().Replace("(+", "(");

		Master *m = gd->master_list.Store(&line.Token());

		m->type = EX;
		m->e  = gd->element_list.Store(&line.Token(), false);
		m->name = line.Token();

		if (!line.NextToken())
			return false;

		RETURN_TYPE tt = line.TokenType();

		if (tt != _upper_ && line.Token().IsAmong(0, "[") && line.Token().Compare("e-", true) != 0)
			return false;

		m->s = gd->species_list.Search(&line.Token());
		if (m->s == NULL)
		{
			m->s = gd->species_list.Store(&line.Token(), false); 
			GetZ(line.Token(), m->s->z);
		}

		m->primary = true;

		if (m->name != "E")
		{
			Element * e = gd->element_list.Store(&m->name, false);
			e->gfw = 0.0;
		}
	}
	return true;
}

bool Database::LoadExchangeSpeciesBlock()
{
	options_list.Clear();
	options_list.AddNew("-no_check");
	options_list.AddNew("-check");
	options_list.AddNew("-gamma");
	options_list.AddNew("-mb");
	options_list.AddNew("-mass_balance");
	options_list.AddNew("-log_k");
	options_list.AddNew("-logk");
	options_list.AddNew("-delta_h");
	options_list.AddNew("-deltah");
	options_list.AddNew("-analytical_expression");
	options_list.AddNew("-a_e");
	options_list.AddNew("-ae");
	options_list.AddNew("-mole_balance");
	options_list.AddNew("-davies");
	options_list.AddNew("-offset");
	options_list.AddNew("-llnl_gamma");
	options_list.AddNew("-add_logk");
	options_list.AddNew("-add_log_k");
	options_list.AddNew("-add_constant");

	s = NULL;

	RETURN_TYPE lt;

	for (;;)
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();

		switch (lt)
		{
		case _error_:
			return false;
		case _keyword_:
			return true;
		case _empty_:
			continue;
		case _option_:
			if (!LoadExchangeSpeciesOption())
				return false;
			break;
		case _none_:
			if (!LoadExchangeSpeciesEquation())
				return false;
			break;
		}
	}
	return true;
}

bool Database::LoadExchangeSpeciesOption()
{
	LDBLE offset;
	
	if (s == NULL)
		return false;

	line.NextToken();

	int option_id;
	options_list.Search(&line.Token(), option_id, true, true);

	switch(option_id)
	{
	case 0: //-no_check
		s->check_equation = false;
		break;

	case 1: //-check
		s->check_equation = true;
		break;

	case 2: //-gamma
		s->exch_gflag = 2; 
		if (!GetValueFromNextToken(s->dha) || !GetValueFromNextToken(s->dhb))
			return false;
		if (s->dha == 0 && s->dhb == 0)
    {
			s->dhb = (LDBLE)99.9;
			s->exch_gflag = 1;
    }
		break;

	case 3:  //-mb
	case 4:  //-mass_balance
	case 12: //-mole_balance
		if(!ReadMassBalance(s))
			return false;
		break;

	case 5: //-log_k
	case 6: //-logk
		if (!GetValueFromNextToken(s->logk[0]))
			return false;
		break;

	case 7: //-delta_h
	case 8: //-deltah
		if (!ReadDeltaH(s->logk[1], s->original_units))
			return false;
		break;

	case 9:  //-analytical_expression
	case 10: //-a_e
	case 11: //-ae
		if (!ReadAnalyticalExpression(&(s->logk[2])))
			return false;
		break;

	case 13: //-davies
    s->exch_gflag = 1;
    s->dha = 0;
    s->dhb = (LDBLE)99.9;
		break;

	case 14: //-offset
		if (!GetValueFromNextToken(offset))
			return false;
		s->logk[0] += offset;
		break;

	case 15: //-llnl_gamma
		s->exch_gflag = 7;		// llnl D-H
		if (!GetValueFromNextToken(s->dha))
			return false;
		break;

	case 16: //-add_logk
	case 17: //-add_log_k
		if(!ReadAddLogK(s->add_logk))
			return false;
		break;

	case 18: //-add_constant
		if (!ReadAddConstant(s->add_logk))
			return false;
		break;
	}

	return true;
}

bool Database::LoadExchangeSpeciesEquation()
{
	int i;
	bool r = true;
	
	s = NULL;
	eos_list->Clear();
	
	if ((r = ParseEq(line.text)))
	{
		rt.token_list->Element(0)->s = gd->species_list.Store(&rt.token_list->Element(0)->name);
		rt.token_list->Element(0)->s->z = rt.token_list->Element(0)->z;
		for (i = 1; i < rt.token_list->Count(); i++)
		{
			rt.token_list->Element(i)->s = gd->species_list.Store(&rt.token_list->Element(i)->name, false);
			rt.token_list->Element(i)->s->z = rt.token_list->Element(i)->z;
		}

		s = rt.token_list->Element(0)->s;

		eos_list->CopyTo(s->eos_list);
		
		ElementOfSpecies *eos;
		for (i = 0; i < eos_list->Count(); i++)
		{
			eos = (*eos_list)[i];
			if (eos->e->name == "C")
				s->carbon = eos->coef;
			else if (eos->e->name == "H")
				s->h = eos->coef;
			else if (eos->e->name == "O")
				s->o = eos->coef;
		}

    LDBLE exchange_coef = 0.0;
    for (i = 1; i < rt.token_list->Count(); i++)
			if (rt.token_list->Element(i)->s->type == EX)
				exchange_coef = rt.token_list->Element(i)->coef;
		s->equiv = exchange_coef;

		rt.CopyReactionsToExchange(s->rxn);

		s->type = EX;
		s->gflag = 4;
    s->exch_gflag = 3;
    s->dha = 0.0;
		s->dhb = 0.0;

    //Save as a phase for inverse modeling only
		Phase *p = gd->phase_list.Store(&s->name);
    if (p != NULL)
		{
			p->formula = s->name;
			p->check_equation = false;
			p->type = EX;

			for (i = 0; i < s->eos_list->Count(); i++)		
				p->eos_list->AddNew((*s->eos_list)[i]); //Will save a COPY
		
			s->rxn->CopyTo(p->rxn);
		}
		else
			r = false;
	}

	return r;
}

bool Database::LoadSurfaceMasterSpeciesBlock()
{
	RETURN_TYPE lt;

	for (;;)
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();

		switch (lt)
		{
		case _error_:
			return false;
		case _keyword_:
			return true;
		case _empty_:
			continue;
		}

		if (!line.NextToken())
			return false;

		line.Token().Replace("(+", "(");
		Master *m = gd->master_list.Store(&line.Token());

		m->type = SURF;
		m->e = gd->element_list.Store(&line.Token(), false);

		if (!line.NextToken() && line.TokenType() != _upper_ && line.Token()[0] != '[')
			return false;

		m->s = gd->species_list.Search(&line.Token());
		if (m->s == NULL) 
		{
			m->s = gd->species_list.Store(&line.Token(), false); //, z, false);
			GetZ(line.Token(), m->s->z);
		}

		m->primary = true;
		String tok = m->name;
		tok.Replace("_", " ");
		tok = tok(0, " ");
		tok += "_psi";
		if (!AddPSIMasterSpecies(tok))
			return false;
	}

	return true;
}

bool Database::LoadSurfaceSpeciesBlock()
{
	options_list.Clear();
	options_list.AddNew("-no_check");
	options_list.AddNew("-check");
	options_list.AddNew("-mb");
	options_list.AddNew("-mass_balance");
	options_list.AddNew("-log_k");
	options_list.AddNew("-logk");
	options_list.AddNew("-delta_h");
	options_list.AddNew("-deltah");
	options_list.AddNew("-analytical_expression");
	options_list.AddNew("-a_e");
	options_list.AddNew("-ae");
	options_list.AddNew("-mole_balance");
	options_list.AddNew("-offset");
	options_list.AddNew("-add_logk");
	options_list.AddNew("-add_log_k");
	options_list.AddNew("-add_constant");
	options_list.AddNew("-cd_music");
	options_list.AddNew("-music");

	s = NULL;

	RETURN_TYPE lt;

	for (;;)
	{
		if (!file->GetLine())
			return false;

		lt = line.LineType();

		switch (lt)
		{
		case _error_:
			return false;
		case _keyword_:
			return true;
		case _empty_:
			continue;
		case _option_:
			if (!LoadSurfaceSpeciesOption())
				return false;
			break;
		case _none_:
			if (!LoadSurfaceSpeciesEquation())
				return false;
			break;
		}
	}
	return true;
}

bool Database::LoadSurfaceSpeciesOption()
{
	if (s == NULL)
		return false;

	if (!line.NextToken())
		return false;

	int option_id;
	options_list.Search(&line.Token(), option_id, true, true);

	LDBLE offset;

	switch(option_id)
	{
	case 0: //-no_check
		s->check_equation = false;
		break;

	case 1: //-check
		s->check_equation = true;
		break;

	case 2:  //-mb
	case 3:  //-mass_balance
	case 11: //-mole_balance
		if(!ReadMassBalance(s))
			return false;
		break;

	case 4: //-log_k
	case 5: //-logk
		if (!GetValueFromNextToken(s->logk[0]))
			return false;
		break;

	case 6: //-delta_h
	case 7: //-deltah
		if (!ReadDeltaH(s->logk[0], s->original_units))
			return false;
		break;

	case 8:  //-analytical_expression
	case 9: //-a_e
	case 10: //-ae
		if (!ReadAnalyticalExpression(&(s->logk[2])))
			return false;
		break;

	case 12: //-offset
		if (!GetValueFromNextToken(offset))
			return false;
		s->logk[0] += offset;
		break;

	case 13: //-add_logk
	case 14: //-add_log_k
		if(!ReadAddLogK(s->add_logk))
			return false;
		break;

	case 15: //-add_constant
		if (!ReadAddConstant(s->add_logk))
			return false;
		break;

	case 16: //-cd_music
	case 17: //-music
		LDBLE cd_music;
    for (int i = 0; i < 5; i++)
		{
			if (!GetValueFromNextToken(cd_music))
				return false;

			s->cd_music[i] = cd_music;
    }

    s->dz[0] = s->cd_music[0] + s->cd_music[3] * s->cd_music[4];
    s->dz[1] = s->cd_music[1] + (1 - s->cd_music[3]) * s->cd_music[4];
    s->dz[2] = s->cd_music[2];
		break;

	default:
		return false;
	}

	return true;
}

bool Database::LoadSurfaceSpeciesEquation()
{
	int i;
	bool r = true;
	
	s = NULL;
	eos_list->Clear();

	if ((r = ParseEq(line.text)))
	{
		rt.token_list->Element(0)->s = gd->species_list.Store(&rt.token_list->Element(0)->name);
		rt.token_list->Element(0)->s->z = rt.token_list->Element(0)->z;
		for (i = 1; i < rt.token_list->Count(); i++)
		{
			rt.token_list->Element(i)->s = gd->species_list.Store(&rt.token_list->Element(i)->name, false);
			rt.token_list->Element(i)->s->z = rt.token_list->Element(i)->z;
		}

		s = rt.token_list->Element(0)->s;
		eos_list->CopyTo(s->eos_list);

		ElementOfSpecies *eos;
		for (i = 0; i < eos_list->Count(); i++)
		{
			eos = (*eos_list)[i];
			if (eos->e->name == "C")
				s->carbon = eos->coef;
			else if (eos->e->name == "H")
				s->h = eos->coef;
			else if (eos->e->name == "O")
				s->o = eos->coef;
		}

    s->equiv = 0.0;
		ReactionToken *r_token;
    for (i = 1; i < rt.token_list->Count(); i++)
		{
			r_token = rt.token_list->Element(i);

			if (r_token->s->type == SURF)
				s->equiv = r_token->coef;
		}

		if (s->equiv == 0.0)
			s->equiv = 1.0;

		Reaction *reaction = s->rxn;
		reaction->token_list->Clear();
		reaction->token_list->SetNewCapacity(rt.token_list->Count());

		long i;
		ReactionToken *rt_ptr, *rt_new;
		for (i = 0; i < rt.token_list->Count(); i++)
		{
			rt_ptr = rt.token_list->Element(i);
			rt_new = reaction->token_list->AddNew();

			rt_new->s = rt_ptr->s;
			rt_new->coef = rt_ptr->coef;
		}

		s->type = SURF;
		s->gflag = 6;
    s->dha = 0.0;
    s->dhb = 0.0;
	}
	else 
		r = false;

	return r;
}

bool Database::TidyDatabase()
{
	gd->species_list.Sort(0, false);
	gd->master_list.Sort(0, false);
	gd->element_list.Sort(0, false);
	gd->phase_list.Sort(0, true);


	//d.PrintMasterList("master_before_tidy_species.txt", gd);
	//d.PrintSpeciesList("species_before_tidy_species.txt", gd);

	if (gd->species_list.Count() > 0 && !TidySpecies())
		return false;

	//d.PrintMasterList("master_after_tidy_species.txt", gd);
	//d.PrintSpeciesList("species_after_tidy_species.txt", gd);

	if (gd->phase_list.Count() > 0 && !TidyPhases())
		return false;

	if (!ComputeGFW(String("H2O"), gfw_water))
		return false;

	gfw_water *= (LDBLE)0.001;

	if (!CheckObrigatorySpecies())
		return false;

	return true;
}

//==============================================================================
// for debug purposes
//==============================================================================

void Database::PrintToFileSpeciesInfo(const char *file_name)
{
	FILE *f = fopen(file_name, "w");

	Species *s;
	for(int i = 0; i < gd->species_list.Count(); i++)
	{
		s = gd->species_list[i];

		fprintf(f, "species name: %s", s->name);
		fprintf(f, " (rxn count: %d)\n", s->rxn->token_list->Count());
	}

	fclose(f);
}

void Database::PrintToFileSpeciesInfo2(const char *file_name)
{
	FILE *saida = fopen(file_name, "w");

	fprintf(saida, "species count: %d\n\n", gd->species_list.Count());
	for (int n = 0; n < gd->species_list.Count(); n++)
	{
		fprintf(saida, "%d: %s\n", n, gd->species_list[n]->name);

		if (gd->species_list[n]->primary)
			fprintf(saida, "   primary: %s\n", gd->species_list[n]->primary->name.CharPtr());
		else
			fprintf(saida, "   primary: <NULL>\n");

		if (gd->species_list[n]->secondary)
			fprintf(saida, "   secondary: %s\n", gd->species_list[n]->secondary->name.CharPtr());
		else
			fprintf(saida, "   secondary: <NULL>\n");

		if (gd->species_list[n]->rxn->token_list->Count() > 0)
		{
			fprintf(saida, "   rxn tokens: ");
			for(int n_t = 0; n_t < gd->species_list[n]->rxn->token_list->Count(); n_t++)
			{
				fprintf(saida, "%s", gd->species_list[n]->rxn->token_list->Element(n_t)->name.CharPtr());
				if (gd->species_list[n]->rxn->token_list->Element(n_t)->s)
					fprintf(saida, "[s: %s]  ", gd->species_list[n]->rxn->token_list->Element(n_t)->s->name.CharPtr());
				else
					fprintf(saida, "[s: null]  ");
			}
			fprintf(saida, "\n");
		}
		else
			fprintf(saida, "   rxn tokens: \n");
	}	

	fclose(saida);
}

void Database::PrintToFilePhasesInfo(const char *file_name)
{
	FILE *saida = fopen(file_name, "w");

	fprintf(saida, "phases count: %d\n\n", gd->phase_list.Count());
	for (int n = 0; n < gd->phase_list.Count(); n++)
	{
		fprintf(saida, "%d: %s\n", n, gd->phase_list[n]->name);

		if (gd->phase_list[n]->rxn->token_list->Count() > 0)
		{
			fprintf(saida, "   rxn tokens: ");
			for(int n_t = 0; n_t < gd->phase_list[n]->rxn->token_list->Count(); n_t++)
			{
				fprintf(saida, "%s", gd->phase_list[n]->rxn->token_list->Element(n_t)->name.CharPtr());
				if (gd->phase_list[n]->rxn->token_list->Element(n_t)->s)
					fprintf(saida, "[s: %s]  ", gd->phase_list[n]->rxn->token_list->Element(n_t)->s->name.CharPtr());
				else
					fprintf(saida, "[s: null]  ");
			}
			fprintf(saida, "\n");
		}
		else
			fprintf(saida, "   rxn tokens: \n");
	}	

	fclose(saida);
}