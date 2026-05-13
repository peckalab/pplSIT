import json
import os
from functools import lru_cache


_ACTIVE_SOUND_ORDER = (
    "background",
    "target",
    "distractor1",
    "distractor2",
    "distractor3",
)
_PASSIVE_PREFIX_ORDER = ("F", "D", "I", "L")
_SOUND_ATTR_KEYS = ("freq", "amp", "duration", "harmonics", "channels", "enabled")


def get_experiment_type(parameters):
    return str(parameters.get("experiment", {}).get("experiment_type", ""))


def detect_session_paradigm(parameters):
    sounds = parameters.get("sound", {}).get("sounds", {})
    return "active" if "background" in sounds else "passive"


def session_cfg_path(src_path, animal, session):
    return os.path.join(src_path, animal, session, f"{session}.json")


@lru_cache(maxsize=None)
def load_session_parameters(src_path, animal, session):
    with open(session_cfg_path(src_path, animal, session), "r") as f:
        return json.load(f)


@lru_cache(maxsize=None)
def get_session_paradigm(src_path, animal, session):
    parameters = load_session_parameters(src_path, animal, session)
    return detect_session_paradigm(parameters)


def _catalog_entry(name, role, spec=None):
    entry = {"name": name, "role": role}
    if isinstance(spec, dict):
        for key in _SOUND_ATTR_KEYS:
            if key in spec:
                entry[key] = spec[key]
    return entry


def _iter_passive_sound_names(sounds):
    excluded = {"noise", "silence"}
    ordered = []

    for prefix in _PASSIVE_PREFIX_ORDER:
        ordered.extend(
            sorted(
                name
                for name, spec in sounds.items()
                if name not in excluded
                and isinstance(spec, dict)
                and spec.get("enabled", True)
                and name.upper().startswith(prefix)
            )
        )

    ordered_set = set(ordered)
    ordered.extend(
        sorted(
            name
            for name, spec in sounds.items()
            if name not in excluded
            and name not in ordered_set
            and isinstance(spec, dict)
            and spec.get("enabled", True)
        )
    )
    return ordered


def build_stimulus_catalog(parameters):
    sounds = parameters.get("sound", {}).get("sounds", {})
    paradigm = detect_session_paradigm(parameters)

    catalog = {
        "-1": _catalog_entry("noise", "noise", sounds.get("noise")),
        "0": _catalog_entry("silence", "silence", sounds.get("silence")),
    }

    if paradigm == "active" or "background" in sounds or "target" in sounds:
        sound_id = 1
        for name in _ACTIVE_SOUND_ORDER:
            spec = sounds.get(name)
            if not isinstance(spec, dict):
                continue
            if name.startswith("distractor") and not spec.get("enabled", False):
                continue

            role = "distractor" if name.startswith("distractor") else name
            catalog[str(sound_id)] = _catalog_entry(name, role, spec)
            sound_id += 1
        return catalog

    for sound_id, name in enumerate(_iter_passive_sound_names(sounds), start=1):
        catalog[str(sound_id)] = _catalog_entry(name, "stimulus", sounds.get(name))

    return catalog
