"""
adapted from Tom Goddard's `bigalpha.py` script
https://rbvi.github.io/chimerax-recipes/big_alphafold/bigalpha.html

open ChimeraX
open this python file
now the command 'bigalpha' will be available

file names should have the format "prefix_F1.pdb, prefix_F2.pdb, ... prefix_FN.pdb"
for N fragments

cd /directory/to/models
bigalpha prefix

"""

def open_multifile_alphafold_model(session, prefix, directory=".",
                                   combine=True, residues_per_file=1400,
                                   overlap=200, align_span=5):
    from os import listdir

    all_filenames = listdir(directory)
    filenames = [
        filename for filename in all_filenames
        if filename.startswith(prefix) and filename.endswith(".pdb")
    ]
    nfiles = len(filenames)
    print("Model is split into %d PDB files." % (nfiles))

    filename = "%s_F%%d.pdb" % prefix

    # find next available model number
    model_id = max([m.id[0] for m in session.models], default=0) + 1

    # open the overlapping component models
    models = []

    fstep = (residues_per_file // overlap) - 1

    for i in range(1, nfiles+1):
        open_next_model(filename % i,
                        residues_per_file, overlap, align_span, models)

    # combine the segment models into one
    model_ids = "#" + ",".join(m.id_string for m in models)

    from chimerax.core.commands import run

    if combine:
        shift_residue_numbers(models)
        copy = run(session, ('combine %s name %s modelId %#d close true'
                             % (model_ids, prefix, model_id)))
        # change the chain IDs for every chain in combined model
        for chain in copy.chains:
            session.logger.info("Mapping back chain %s to A" % (chain.chain_id))
            chain.chain_id = "A"

        models = [copy]

    else:
        run(session, "rename %s %s id #%d" % (model_ids, prefix, model_id))

    run(session, "light full")
    run(session, "view")

    return models


def open_next_model(path, residues_per_file, overlap, align_span, models):
    from chimerax.core.commands import run
    model = run(session, 'open %s' % path)[0]
    id = model.id_string
    run(session, 'color bfactor #%s palette alphafold' % id)
    run(session, 'hide #%s cartoon ; show #%s atoms ; style #%s sphere' % (id,id,id))

    end_match = '%d-%d' % (residues_per_file-5, residues_per_file-1)  # 1395-1399
    if models:
        last_model = models[-1]
        run(session, 'align #%s:%d-%d@CA to #%s:%s@CA'
            % (model.id_string, overlap-align_span+1, overlap,
               last_model.id_string, end_match))
        run(session, 'delete #%s:1-%d' % (model.id_string, overlap))

        # form peptide bond between  models
        # #<previous_model>:<overlap>@C
        # #<last_model>:<overlap+1>@N
        session.logger.info("#############  DEBUG ###############")
        session.logger.info(
            "%s:%d@C <---> %s:%d@N" % (
                last_model.id_string, residues_per_file, model.id_string, overlap+1)
        )

        run(session, "build join peptide #%s:%d@C #%s:%d@N" % (
            last_model.id_string, residues_per_file, model.id_string, overlap+1
        ))

    models.append(model)


def shift_residue_numbers(structures):
    # First adjust residue numbers of each segment
    rnext = 1
    for i,m in enumerate(structures):
        rnums = m.residues.numbers
        rmax, rmin = rnums.max(), rnums.min()
        m.residues.numbers += rnext - rmin
        rnext += rmax-rmin+1


def register_command(session):
    from chimerax.core.commands import CmdDesc, register, StringArg, OpenFolderNameArg, BoolArg
    desc = CmdDesc(required=[('prefix', StringArg)],
                   keyword=[('directory', OpenFolderNameArg),
                            ('combine', BoolArg)],
                   synopsis='Open multifile AlphaFold model')
    register('bigalpha', desc, open_multifile_alphafold_model, logger=session.logger)

register_command(session)
