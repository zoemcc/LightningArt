using DrWatson
@quickactivate "LightningArt"
using LightningArt
using Luxor
using LinearAlgebra
using ForwardDiff
import GeometryBasics
GB = GeometryBasics
using LightGraphs, MetaGraphs
using FileIO
using Revise
using ImageMagick
using ColorTypes, ColorVectorSpace
using Plots, GraphPlot
using Interpolations

guitar_img = load(datadir("exp_raw/guitar_model_1.png"))
plot(guitar_img)
bool_mask_to_blackwhite(bool) = bool .* RGBA(0.0, 0.0, 0.0, 1.0) .+ !bool .* RGBA(1.0, 1.0, 1.0, 1.0)

binary_collision_mask = map(col->norm(col) < 2, guitar_img)

function generate_edge_graph(binary_collision_mask)

    height, width = size(binary_collision_mask)
    pixel_graph = LightGraphs.SimpleGraph(height * width)
    linear_indices_mask = LinearIndices(binary_collision_mask)
    for j in 1:width-1, i in 1:height-1
        if !binary_collision_mask[i, j] && !binary_collision_mask[i+1, j]
            add_edge!(pixel_graph, linear_indices_mask[i, j], linear_indices_mask[i+1, j])
        end

        if !binary_collision_mask[i, j] && !binary_collision_mask[i, j+1]
            add_edge!(pixel_graph, linear_indices_mask[i, j], linear_indices_mask[i, j+1])
        end
    end
    pixel_graph
end

pixel_graph = generate_edge_graph(binary_collision_mask)

connected_components_guitar = LightGraphs.connected_components(pixel_graph)
large_connected_components_guitar = filter(x-> length(x) > 1, connected_components_guitar)
linear_indices_mask = LinearIndices(binary_collision_mask)
index_of_an_interior_pixel = linear_indices_mask[200, 100]
inside_connected_guitar_indices = Set(filter(x->index_of_an_interior_pixel in x, large_connected_components_guitar)[1])
length(inside_connected_guitar_indices)
inside_guitar_mask = map(linearindex->linearindex in inside_connected_guitar_indices, linear_indices_mask)
plot(bool_mask_to_blackwhite.(inside_guitar_mask))

inside_guitar_interpolation = interpolate(inside_guitar_mask, BSpline(Constant()))
function inside_guitar_function(p) 
    interp_coords = clamp.(p, 0.0, 1.0) .* (size(inside_guitar_mask) .- 1) .+ 1
    inside_guitar_interpolation(interp_coords...) ≈ 1
end

height, width = size(inside_guitar_mask)
points_box = [GB.Point2(x, y) for x in 0.0:1/height:1.0, y in 0.0:1/width:1.0]
guitar_interp_mask = inside_guitar_function.(points_box)
plot(bool_mask_to_blackwhite.(guitar_interp_mask))




plot(bool_mask_to_blackwhite.(binary_collision_mask))


### TODO: define a function that is 1 inside the guitar and 0 outside