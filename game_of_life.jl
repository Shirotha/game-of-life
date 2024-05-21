include("data.jl")
include("display.jl")

function render(renderer::MatrixRenderer)
    display(renderer)
    print("\e[$(size(renderer.data, 1)+1)A")
    sleep(1/10)
end

function main(width, height)
    a = Domain{Bool}(width, height; stencil_size = 3)
    b = similar(a)
    boundary = (Periodic(), Periodic())

    renderer = BoolRenderer(view(a))

    update!(a) do s
        rand(Bool)
    end
    bound!(a, boundary...)

    render(renderer)

    Base.exit_on_sigint(false)
    try
    while true
        update!(a, b) do s
            neighbours = count(s) - s[2, 2]
            if s[2, 2]
                return 2 <= neighbours <= 3
            else
                return neighbours == 3
            end
        end
        bound!(b, boundary...)

        a, b = b, a

        render(renderer)
    end
    catch ex
        if !(ex isa InterruptException)
            show(ex)
        end
    finally
        Base.exit_on_sigint(true)
    end
end

main(10, 10)
