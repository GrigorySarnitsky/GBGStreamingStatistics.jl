# HDF5 persistence for the summaries — loaded automatically when both GBGStreamingStatistics and HDF5 are
# present (a weak-dependency package extension, so the core stays free of the HDF5 binary). A summary maps
# to a self-describing HDF5 *group* (arrays -> datasets, scalars -> attributes, plus a type/version tag),
# never a single dataset; the consumer owns the file, so every method takes an open `File`/`Group`.
module GBGStreamingStatisticsHDF5Ext

using HDF5
using GBGStreamingStatistics: Histogram1D, Histogram2D
import GBGStreamingStatistics.Histograms: write_histogram, read_histogram

const TAG_1D = "GBGStreamingStatistics.Histogram1D"
const TAG_2D = "GBGStreamingStatistics.Histogram2D"
const SCHEMA_VERSION = 1


# --- write: create the subgroup, store the host arrays as datasets and the scalars as attributes ---

function write_histogram(parent::HDF5.H5DataStore, name::AbstractString, h::Histogram1D)
  g = create_group(parent, name)
  g["counts"] = Array(h.counts)              # Array(...) copies a device array to the host; no-op on CPU
  g["powersums"] = Array(h.powersums)
  a = attrs(g)
  a["type"] = TAG_1D
  a["version"] = SCHEMA_VERSION
  a["lower"] = h.lower
  a["upper"] = h.upper
  a["nbins"] = length(h.counts)              # stored for self-description; reconstruction uses the array
  a["order"] = size(h.powersums, 1)
  return g
end


function write_histogram(parent::HDF5.H5DataStore, name::AbstractString, h::Histogram2D)
  g = create_group(parent, name)
  g["counts"] = Array(h.counts)
  g["powersums"] = Array(h.powersums)
  a = attrs(g)
  a["type"] = TAG_2D
  a["version"] = SCHEMA_VERSION
  a["xlower"] = h.lower_x
  a["xupper"] = h.upper_x
  a["nx"] = h.n_x
  a["ylower"] = h.lower_y
  a["yupper"] = h.upper_y
  a["ny"] = h.n_y
  a["lag"] = h.lag
  return g
end


# the sugar: `parent[name] = summary` mirrors HDF5.jl's own dataset-writing `setindex!`. Dispatch on the
# concrete `File`/`Group` (HDF5's own parent type) and the specific value type, so this method is strictly
# more specific than — never ambiguous with — HDF5's array-writing `setindex!(::Union{File,Group}, val, path)`.
function Base.setindex!(parent::Union{HDF5.File,HDF5.Group}, h::Union{Histogram1D,Histogram2D}, name::AbstractString)
  write_histogram(parent, name, h)
  return parent
end


# --- read: dispatch on the stored type tag, reconstruct via the public constructor (which fixes the
# element types from the stored datasets, so the count-eltype cone variant is preserved), copy arrays in ---

function read_histogram(parent::HDF5.H5DataStore, name::AbstractString)
  haskey(parent, name) ||
    throw(ArgumentError("read_histogram: no object named \"$name\" under the given parent"))
  g = parent[name]
  g isa HDF5.Group ||
    throw(ArgumentError("read_histogram: \"$name\" is a dataset, not a histogram group"))
  haskey(attrs(g), "type") || throw(ArgumentError(
    "read_histogram: group \"$name\" has no `type` attribute — it was not written by write_histogram"))
  tag = attrs(g)["type"]
  if tag == TAG_1D
    counts = read(g["counts"])
    powersums = read(g["powersums"])
    h = Histogram1D(attrs(g)["lower"], attrs(g)["upper"], length(counts);
                    order = size(powersums, 1), accumulator_type = eltype(powersums), count_type = eltype(counts))
    copyto!(h.counts, counts)
    copyto!(h.powersums, powersums)
    return h
  elseif tag == TAG_2D
    counts = read(g["counts"])
    powersums = read(g["powersums"])
    h = Histogram2D(attrs(g)["xlower"], attrs(g)["xupper"], Int(attrs(g)["nx"]),
                    attrs(g)["ylower"], attrs(g)["yupper"], Int(attrs(g)["ny"]);
                    lag = Int(attrs(g)["lag"]), accumulator_type = eltype(powersums), count_type = eltype(counts))
    copyto!(h.counts, counts)
    copyto!(h.powersums, powersums)
    return h
  else
    throw(ArgumentError("read_histogram: unrecognized histogram type tag \"$tag\" in group \"$name\""))
  end
end

end # module GBGStreamingStatisticsHDF5Ext
