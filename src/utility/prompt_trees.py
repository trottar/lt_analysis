#! /usr/bin/python

def normalize_epsset(epsset):
    eps_text = str(epsset).strip().lower()
    if eps_text.endswith("e"):
        eps_text = eps_text[:-1]
    return eps_text


def use_prompt_rf_tree(epsset):
    return normalize_epsset(epsset) == "low"


def get_prompt_tree_name(particle_type, epsset):
    rf_suffix = "RF" if use_prompt_rf_tree(epsset) else "noRF"
    return "Cut_{}_Events_prompt_{}".format(particle_type.capitalize(), rf_suffix)


def get_subtraction_prompt_tree_name(particle_type, epsset):
    particle_text = str(particle_type).strip().lower()
    if particle_text == "pion":
        rf_suffix = "noRF"
    else:
        rf_suffix = "RF" if use_prompt_rf_tree(epsset) else "noRF"
    return "Cut_{}_Events_prompt_{}".format(particle_type.capitalize(), rf_suffix)
