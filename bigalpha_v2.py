def open_multifile_alphafold_model(
        session, prefix, directory=".",
        combine=True,
        residues_per_file=1400,
        overlap=200,
        align_span=5,
):
    from os import listdir, path
    import re

    all_filenames = listdir(directory)

    # filename pattern to look for
    fptn = r"^%s_seg(\d+)_.*\.pdb$" % (prefix)
    ptn = re.compile(fptn)

    segment_order = []
    pdb_files = []

    for f in all_filenames:
        m = ptn.match(f)
        if m:
            segment_order.append(int(m.groups()[0]))
            pdb_files.append(f)

    # sort pdb files based on segment number
    sorted_pdb_files = [f for _, f in sorted(zip(segment_order, pdb_files))]
    nfiles = len(sorted_pdb_files)

    # find next available model number
    model_id = max([m.id[0] for m in session.models], default=0) + 1

    # open the overlapping segments
    models = []

    for f in sorted_pdb_files:
        full_fn = path.join(directory, f)
        open_next_model(
            full_fn, residues_per_file, overlap, align_span, models
        )

    from chimerax.core.commands import run

    if combine:
        model_ids = "#" + ",".join(m.id_string for m in models)

        for m in models:
            m.chain_id = "A"

        # before we combine models, get the residues to join the termini
        c_terms = []
        n_terms = []

        for n in range(len(models)-1):
            c_terms.append(models[n].residues[-1].number)
            n_terms.append(models[n+1].residues[0].number)

        # combine models by joining their termini
        model_numbers = sorted([int(m.id_string) for m in models])
        root_model = model_numbers[0]

        for i, (c_term, n_term) in enumerate(
                zip(c_terms, n_terms), start=model_numbers[1]
        ):
            run(session,
                "build join peptide #%d:%d@C #%d:%d@N" % (
                    root_model,
                    c_term,
                    i,
                    n_term
                ))

        run(session, "rename #%s %s" % (str(root_model), prefix))
        run(session, "light soft")
        run(session, "view")


def open_next_model(fn, residues_per_file, overlap, align_span, models):
    from chimerax.core.commands import run

    model = run(session, "open %s" % fn)[0]
    model_id = model.id_string

    # adjust appearance
    run(session, "color bfactor #%s palette alphafold" % model_id)
    run(session,
        "hide #%s cartoon; show #%s atoms; style #%s sphere" %
        (model_id, model_id, model_id))

    # if the first segment has already been loaded:
    if models:
        last_model = models[-1]
        last_residue = last_model.residues[-1].number

        # align current model to the last one by 'align_span'sel 
        run(
            session,
            "align #%s:%d-%d@CA to #%s:%d-%d@CA" % (
                model.id_string,
                overlap - align_span + 1,
                overlap,
                last_model.id_string,
                last_residue - align_span + 1,
                last_residue
            )
        )

        # replace bfactors of last model with the current segment
        for i in range(-align_span, 0):
            new_bfactor = model.residues[overlap + i].atoms[0].bfactor
            session.logger.info("changing bfactor for residue %d from %f to %f" % (
                models[-1].residues[i].number,
                models[-1].residues[i].atoms[0].bfactor,
                new_bfactor
            ))
            for atom in models[-1].residues[i].atoms:
                atom.bfactor = new_bfactor


        # and delete the overlapping part of current model
        run(
            session,
            "delete #%s:1-%d" % (model.id_string, overlap)
        )

        # renumber the residue of current model (after deletion)
        model.renumber_residues(model.residues, last_residue + 1)

        # join the C-term of previous model to N-term of current model
        # run(
        #     session,
        #     "build join peptide #%s:%d@C #%s:%d@N" % (
        #         last_model.id_string,
        #         residues_per_file,
        #         model.id_string,
        #         overlap + 1
        #     )
        # )

    models.append(model)


def register_command(session):
    from chimerax.core.commands import (
        CmdDesc, register, StringArg, OpenFolderNameArg, BoolArg, IntArg
    )
    desc = CmdDesc(required=[('prefix', StringArg)],
                   keyword=[('directory', OpenFolderNameArg),
                            ('combine', BoolArg),
                            ('residues_per_file', IntArg),
                            ('overlap', IntArg),
                            ('align_span', IntArg)],
                   synopsis='Open multifile AlphaFold model')
    register('bigalpha', desc, open_multifile_alphafold_model,
             logger=session.logger)


register_command(session)
