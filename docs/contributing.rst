.. _contributing:


Contributing to the GDT
=======================
Community contributions to the GDT are welcome and make the GDT a better 
toolkit for everyone to use.  There are various types of contributions.  These 
include:

* Bug fixes
* New features
* Documentation improvements
* API improvements
* New mission packages

For each contribution, we consider it best practice to first open an issue (with
an appropriate ``label``).  This may be as simple as notifying us (and other 
users!) that there is a bug.  Perhaps you also have a proposed solution to 
fixing the issue.  If so, then you can outline your proposed fix when you open 
the issue, and we will give you feedback about whether we think that would be a
useful and appropriate fix.

**Simple typographical fixes and clarifications made in the documentation do not require
the creation of an issue.**

If your proposed modifications are significant (e.g. sizable document change, API improvements,
new feature), we highly recommend that you detail the propose change in the issue and wait for feedback.
This is to save you precious time in the event that we decide not to accept your proposed solution.
Often the proposed solution can break functionality elsewhere or can be simplified, and we would like to have
a chance to provide useful feedback before you begin coding. Waiting for feedback
is not a requirement, but merely reduces the chance of additional changes prior to having
your pull request accepted.

If you are submitting code modifications, we require that you create a unit test to confirm
expected operation of the code if those modifications aren't already covered by an
existing unit test.

The usual sequence of events are:

1. Create an issue describing the proposed changes.
2. Waiting for feedback if desired.
3. Fork from main branch.
4. Use your fork to add your changes to the code base.
5. Create unit tests showing that your changes work as intended.
6. Create `Pull Request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_ with a comment explaining how it closes the issue you created.


Opening an Issue
-----------------
As mentioned above, all types of issues are welcome; however, an issue without
sufficient description is not helpful.  When opening an issue, please provide
detail about the issue.  If the issue is a bug, include a description of the bug,
how it was encountered, and a stack trace is often helpful.  If there is a 
proposed fix, explain what that fix is and where it would be implemented.  For
feature requests and API improvements, detail needs to be provided about what 
the new feature/improvement is and how it might be implemented.


Requisites for Pull Requests
----------------------------
For all pull requests (PRs), there needs to be sufficient detail of the 
issue/feature and the corresponding implementation, including rationale where
appropriate. If an issue is opened first, as suggested above, the PR can simply
reference the open issue for these details if the opened issue already contains
a sufficient description.

For bug fixes, a unit test may be required to ensure that the implementation 
does, in fact, fix the bug.  For PRs covering new features and API improvements,
unit tests covering the new features/improvements are *required*, as are 
documentation updates describing the new features or API improvements.  Unit 
tests are important for providing at least a modicum of assurance that the code 
performs as intended, and documentation is necessary because features will 
go unused if no one knows how to use them.  There may be some exceptions to unit
tests, such as code that makes figures.  In which case, the PR should include
attachments showing figures produced by the proposed code.

It is strongly encouraged to run all unit tests that are relevant to your PR 
before submitting it.

Contributing to ``gdt-core``
----------------------------
When contributing to GDT, one must first think about the most appropriate place
for the contribution.  ``gdt-core`` is where all the core functionality and the
base classes that most of the mission packages rely on.  If there is a bug fix
or feature that is related to that functionality (e.g binning, background 
fitting, base data classes, etc.) then that belongs in ``gdt-core``.

Particular care is taken in reviewing fixes/features for ``gdt-core`` since any
change may have impacts to some or all of the missions within the GDT family.
This means that when we consider issues like API improvements, we may balance
this with the impact it would have on all of the mission packages.  Major API
changes, particularly in ``gdt-core``, rarely get accepted for the current release.
Instead it will have to wait for the next planned major version release, as
it could certainly imply changes that propagate to the mission packages.

Contributing to mission packages
--------------------------------
Generally, contributing to a mission package is a little simpler than
contributing to ``gdt-core``.

Contributions to a mission package of a fix or
feature applies to that mission, and that mission alone.  This usually
impacts mission-specific data or functionality that is not part of ``gdt-core``.
In some cases, it may be possible to implement a new feature in a mission 
package and then work to generalize that feature so that it can be implemented
in ``gdt-core`` and available for other missions to use.

Contributing new missions packages
----------------------------------
Users are welcome to create their own mission packages within the GDT namespace
and host on PyPI.  We would be happy to include a reference to the mission 
package in our :ref:`list of supported missions<mission-packages>`, making it clear who is providing
support for that package.  In some cases, we would be willing to host and support
new mission packages under certain conditions.  Those who are interested in that
option should contact STI-Software@usra.edu.

