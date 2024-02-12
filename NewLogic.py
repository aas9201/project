

class Part:
    def __init__(self, component_type, role):
        self.component_type = component_type
        self.role = role

    def form_cargo(self,*args):
        self.cargo_parts = list(args)


promoter = Part("cargo", "promoter")
rbs = Part("cargo", "rbs")
cds = Part("cargo", "cds")

cargo = Part("cargo", "cargo")
cargo.form_cargo(promoter,rbs,cds)
print(cargo.cargo_parts)
for i in cargo.cargo_parts:
    print(i.component_type)

alpha = ['A', 'B', 'F']
beta = ['C', 'D' , 'E']
insert = alpha.index('B') + 1

alpha[insert:insert] = beta
print(alpha)

A=1
B=2
C= None
D=6
E=7
F=8
alpha = [A,B,C,D,E,F]
alpha[2:3] = [3,4,5]
print(alpha)

