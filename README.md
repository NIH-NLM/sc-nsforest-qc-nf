# nsforest-nf

## Docker

See: [nsforest-docker](https://github.com/NIH-NLM/nsforest-docker)

## Singularity

### To build a singularity container

`$ singularity build nsforest_latest.sif docker://ralatsdio/nsforest:latest`

### To push a singularity container

`$ singularity push -U nsforest_latest.sif library://ralatsdio/nlm-kb/nsforest_latest.sif`

### To run a singularity container

`$ singularity shell nsforest_latest.sif`
