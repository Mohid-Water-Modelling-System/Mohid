# Domain Decomposition

The Domaim Decomposition project aim is to enable Mohid's paralelization using MPI. Domain Decomposition divides a domain in a arbitrary number of subdomains that communicate among them. In the frontier between two subdomains border condtions are exchanged reciprocally between them.

## Subdomains

The domain decomposition project will start with the following subdomain charactistics:

1. Subdomains are derived from a domain dividing it through straight lines in the indexes' directions (usually horizontal or vertical lines). Submodels inherit indexes from the initial domain.
