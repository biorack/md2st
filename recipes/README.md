# The Recipe Folder

This folder contains fragments of code that are examples of how compounds had to
be processed at a given time.  Unfortunately, we could not devise an all
encompassing recipe that would handle all cases.  As a temporary solution, we
will tag compound-database entries with code snippets found here.  This folder
is a simple way to associate what we did to structures as given by various
compound structure providers.

# How to use this Folder

When you add new molecules to our compound database, take the recipe you used
and make a new file and add it to this folder.  It does not have to be perfect
the idea is to catch 99% of what you had to do to process the molecules.  The
alternative is having to store all this information in your mind.

If you are not using our compound database then this content will not be very
useful.  

The general idea of these recipes is:

1. Import compounds.  This is often requiring very custom code.  Structures are
not typically stored appropriately for 100% of the compounds we import from a
third party database.
2. Desalt and neutralize.  This is typically following a standard recipe that
can be found on RDKit help pages all over the internet.  Sometimes, this recipe
does not work.  Changes to this recipe need to be captured.  Where possible,
we tried to follow the most basic recipe that is widely posted.
3. Canonical tautomer.  This is the single most important part.  There does is
not a publicly available tool for doing this that works 100% of the time.
Typically, we use the MolVS package.  At the time of writing this, it is the
best tool available.  There are some corner cases with phosphate that are not
right; also if there are stereo-centers on a ring your stereochemistry will be
removed.  If you can tolerate these two problems then MolVS is very useful and
powerful.
4. Look towards the future.  We fully expect that reliable 3d relaxed structures
will be available in the very near term. Modeling structures this way will
solve almost all standardization problems.  There are already algorithms for
comparing 3d structures that work very well. We do not do this right now, but
we want to setup the databases so they can be reprocessed and consolidated
using 3d relaxed structures without losing compound identifiers that are
associated with other observations.

Ben Bowen
May 12, 2020
