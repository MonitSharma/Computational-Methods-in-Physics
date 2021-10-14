from manim import *
class SquareToCircle(Scene):
    def construct(self):
        circle = Circle()
        square = Square()


        self.play(Create(square))
        self.play(Transform(square, circle))
        self.play(FadeOut(square))




class displayText(Scene):
    def construct(self):
        # Create Text objects
        first_line = Text('Partial Differential Equations')
        second_line = Text('using Manim')
        third_line = Text('Computational Physics', color=RED)

        # Position second line underneath first line
        second_line.next_to(first_line, DOWN)

        # Displaying text
        self.wait(1)
        self.play(Write(first_line), Write(second_line))
        self.wait(1)
        self.play(ReplacementTransform(first_line, third_line), FadeOut(second_line))
        self.wait(2)


class displayEquations(Scene):
    def construct(self):
        # Create Tex objects
        first_line = Text('Partial Differential Equation')
        second_line = Text('Wave Equation : Guitar String')
        equation = Tex('$d\\left(p, q\\right)=d\\left(q, p\\right)=\\sqrt{(q_1-p_1)^2+(q_2-p_2)^2+...+(q_n-p_n)^2}=\\sqrt{\\sum_{i=1}^n\\left(q_i-p_i\\right)^2}$')
        
        # Position second line underneath first line
        second_line.next_to(first_line, DOWN)

        # Displaying text and equation
        self.play(Write(first_line), Write(second_line))
        self.wait(1)
        self.play(ReplacementTransform(first_line, equation), FadeOut(second_line))
        self.wait(3)






class CreateGraph(Scene):
    def construct(self):
        axes = Axes(
            x_range=[-3, 3],
            y_range=[-5, 5],
            axis_config={"color": BLUE},
        )

        # Create Graph
        graph = axes.get_graph(lambda x: x**2, color=WHITE)
        graph_label = axes.get_graph_label(graph, label='x^{2}')

        graph2 = axes.get_graph(lambda x: x**3, color=WHITE)
        graph_label2 = axes.get_graph_label(graph2, label='x^{3}')

        # Display graph
        self.play(Create(axes), Create(graph), Write(graph_label))
        self.wait(1)
        self.play(Transform(graph, graph2), Transform(graph_label, graph_label2))
        self.wait(1)


        

class threeDGraph(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        circle=Circle()
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        text3d = Text("This is a 3D text")
        self.add(circle,axes)
        self.add_fixed_in_frame_mobjects(text3d)
        text3d.to_corner(UL)
        self.add(axes)
        self.wait()