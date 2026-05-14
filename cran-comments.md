## Resubmission
This is a resubmission. In this version, I have fully addressed the reviewer's (Konstanze Lauseker) comments:

* Added core literature references with DOIs in the `Description` field of the DESCRIPTION file.
* Replaced `\dontrun{}` with `\donttest{}` for examples requiring > 5s execution time, and fully unwrapped fast-running examples to allow automatic testing.
* Ensured no functions or examples write to the user's home filespace by default. All file-writing examples and functions now strictly use `tempdir()`.

## Test environments
* local Windows 10 Professional (22H2), R 4.5.1
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Reverse dependencies
This is a new release, so there are no reverse dependencies.
