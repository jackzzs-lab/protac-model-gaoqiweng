import click
from schrodinger.structure import StructureReader, StructureWriter

def reader(f_path: str, index: int = 1):
    with StructureReader(f_path, index=index) as reader:
        for st in reader:
            yield st


def read_one(f_path: str, index: int = 1):
    try:
        return next(reader(f_path, index=index))
    except StopIteration:
        return None

def writer(f_path: str):
    return StructureWriter(f_path)


def write_one(st, f_path: str):
    with writer(f_path) as f:
        f.append(st)

@click.command()
@click.argument("inputs", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option("--output", "-o", required=True)
def cli(inputs, output):
    if len(inputs) < 2:
        raise ValueError('two inputs are required.')
    result = None
    for i in inputs:
        if not result:
            result = read_one(i)
        else:
            result = result.merge(read_one(i), copy_props=True)
    write_one(result, output)

if __name__ == "__main__":
    cli()