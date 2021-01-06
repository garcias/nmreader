# NMReader

Python scripts to read and process NMR FID from NMReady dx file.

Module `nmreader` provides a `Spectrum` class that parses the .dx file to reconstruct the FID, then computes the FT spectrum with both frequency and chemical shift axes. Also provides functions to work with intervals (e.g. parsing from a string, baseline correction with interval-based mask).

For use in Jupyter/Colab notebook to provide GUI with ipywidgets.

## Intallation

In a Jupyter/Colab notebook type the following into a cell and run.

```bash
    ! python -m pip install git+https://github.com/garcias/nmreader.git &>> install.log
    ! wget https://raw.githubusercontent.com/garcias/nmreader/main/test.dx
```

Then import `nmreader` and create a `Spectrum` from the test file.

```python
    import nmreader
    file = 'test.dx'
    spec = nmreader.Spectrum( file )
```

Then use your favorite chart package to plot `spec.fft.real` vs `spec.shift`, to see the spectrum on a chemical shift axis.
If you want to see the FID, plot `spec.fid.real` vs `spec.time`.

To try out the phase correction:

```python
    phased = spec.phase( 1490, -10, -25 )  # pivot the first-order correction at index 1490
```

And then plot `phased.real` vs `spec.shift`.
