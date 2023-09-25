This directory holds dict-serialized MolecularSystems,
    which allows referencing FCI and FES energies without needing to re-diagonalize.

At present, this directory needs to be included in the job folder to avoid errors.

Ultimately, the whole MolecularSystems thing should be standardized somehow but meh.

Most frustratingly, Dict is not as stable as I'd have liked:
    older versions of Julia (ie. the one on ARC) can't read objects serialized with newer versions.
Thus, ARC has to fill up the sys folder itself.
Thus, .gitignore should tend to ignore anything in the sys directory,
    lest we put an incompatible serial on ARC.

But we need to still *have* the folder on ARC to avoid another tedious error.
So this README serves as an object I can force-add such that the folder is included in git.