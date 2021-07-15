Hello world

Dracula - A Clinical Trials Drug - Adverse Event Mapper
=====

**Dracula is a tool to be used in combination with the [CTTI AACT database](https://aact.ctti-clinicaltrials.org/). It
does two things: First, it generates an additional table called 'result_group_ingredient' that links the reported adverse events
from the 'result_groups' table to an RxNorm ID of the drugs provided to the specific group for which this adverse event
occurred. Second it adds an additional pt_code column to the reported_events table containing the MedDRA preferred term
code for the specific adverse event.**

Note:
The master branch of this project should work, but the project as whole is still under construction. Although a decent
effort is made to make this mapping decent, Dracula is not perfect. Contributions, comments and suggestions are very
much welcome.
