#! /usr/bin/python

def normalize_epsset(epsset):
    eps_text = str(epsset).strip().lower()
    if eps_text.endswith("e"):
        eps_text = eps_text[:-1]
    return eps_text


def use_prompt_rf_tree(epsset):
    return normalize_epsset(epsset) == "low"


def normalize_rf_state(rf_state=None, epsset=None):
    if rf_state is None:
        return "RF" if use_prompt_rf_tree(epsset) else "noRF"

    state_text = str(rf_state).strip().lower()
    if state_text in ("rf", "prompt_rf", "cut_rf"):
        return "RF"
    if state_text in ("norf", "no_rf", "prompt_norf", "cut_norf"):
        return "noRF"
    raise ValueError(
        "Unsupported RF state '{}'. Expected 'RF' or 'noRF'.".format(rf_state)
    )


def get_prompt_tree_name(particle_type, epsset=None, rf_state=None):
    rf_suffix = normalize_rf_state(rf_state=rf_state, epsset=epsset)
    return "Cut_{}_Events_prompt_{}".format(particle_type.capitalize(), rf_suffix)


def get_rand_tree_name(particle_type, epsset=None, rf_state=None):
    rf_suffix = normalize_rf_state(rf_state=rf_state, epsset=epsset)
    return "Cut_{}_Events_rand_{}".format(particle_type.capitalize(), rf_suffix)
