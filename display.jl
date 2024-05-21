struct MatrixRenderer{T, A <: AbstractMatrix{T}, F}
    data::A
    render_element::F
end

function Base.show(io::IO, ::MIME"text/plain", renderer::MatrixRenderer)
    if get(io, :compact, false)
        show(io, renderer)
    else
        buffer = string()
        linebreak = "\n"
        ys, xs = axes(renderer.data)
        @inbounds for y in ys
            for x in xs
                buffer *= renderer.render_element(renderer.data[y, x])
            end
            buffer *= linebreak
        end
        print(io, buffer)
    end
end

render_bool(x::Bool) = x ? "⬜" : "⬛"
const BoolRenderer{A} = MatrixRenderer{Bool, A, typeof(render_bool)}
BoolRenderer(data) = MatrixRenderer(data, render_bool)

render_binary(threshold) = x -> render_bool(x >= threshold)
const BinaryRenderer = MatrixRenderer
BinaryRenderer(data, threshold) = MatrixRenderer(data, render_binary(threshold))