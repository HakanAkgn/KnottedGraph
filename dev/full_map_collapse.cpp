// Dense-storage implementation of the existing elementary cubical collapses.
// No floating geometry, thinning, smoothing, or topological fallback occurs.
// Queue order agrees with cubical_spine.py: higher parent dimension, then face.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
using Id = uint32_t;
struct Grid {
    std::array<int,3> n;
    std::array<uint64_t,3> stride;
    uint64_t size;
    Grid(int x,int y,int z):n{2*x+1,2*y+1,2*z+1},stride{uint64_t(2*y+1)*(2*z+1),uint64_t(2*z+1),1} {
        size=uint64_t(n[0])*n[1]*n[2];
        if(size>=std::numeric_limits<Id>::max()) throw std::runtime_error("cell index capacity exceeded");
    }
    std::array<int,3> xyz(Id id) const {
        return {int(id/stride[0]),int((id/stride[1])%n[1]),int(id%n[2])};
    }
    Id id(int x,int y,int z) const {return Id(x*stride[0]+y*stride[1]+z);}
    int dim(Id id) const {auto p=xyz(id);return (p[0]&1)+(p[1]&1)+(p[2]&1);}
};
using Cells=std::vector<uint8_t>;
Cells build(const Grid& g,const std::vector<uint8_t>& mask,int nx,int ny,int nz) {
    Cells result(g.size,0);
    for(int x=0;x<nx;++x) for(int y=0;y<ny;++y) for(int z=0;z<nz;++z) {
        if(!mask[(uint64_t(x)*ny+y)*nz+z]) continue;
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) for(int k=0;k<3;++k)
            result[g.id(2*x+i,2*y+j,2*z+k)]=1;
    }
    return result;
}
// Generator predicate. The verifier below deliberately does not call it.
int64_t free_parent(const Grid& g,const Cells& c,Id face) {
    if(!c[face]) return -1;
    auto p=g.xyz(face);int count=0;Id parent=0;
    for(int a=0;a<3;++a) if(!(p[a]&1)) {
        if(p[a]>0 && c[face-g.stride[a]]) {parent=Id(face-g.stride[a]);++count;}
        if(p[a]+1<g.n[a] && c[face+g.stride[a]]) {parent=Id(face+g.stride[a]);++count;}
    }
    if(count!=1) return -1;
    p=g.xyz(parent);
    for(int a=0;a<3;++a) if(!(p[a]&1)) {
        if(p[a]>0 && c[parent-g.stride[a]]) return -1;
        if(p[a]+1<g.n[a] && c[parent+g.stride[a]]) return -1;
    }
    return parent;
}
std::array<uint64_t,4> counts(const Grid& g,const Cells& cells) {
    std::array<uint64_t,4> r{0,0,0,0};
    for(uint64_t i=0;i<g.size;++i) if(cells[i]) ++r[g.dim(Id(i))];
    return r;
}
void replay(const Grid& g,Cells& c,const std::vector<std::array<Id,2>>& moves) {
    for(size_t index=0;index<moves.size();++index) {
        Id f=moves[index][0],p=moves[index][1];
        if(f>=g.size || p>=g.size || !c[f] || !c[p]) throw std::runtime_error("replay missing cell");
        auto fp=g.xyz(f),pp=g.xyz(p);
        int changed=0;
        for(int a=0;a<3;++a) if(fp[a]!=pp[a]) {
            if(!(pp[a]&1) || (fp[a]&1) || std::abs(fp[a]-pp[a])!=1)
                throw std::runtime_error("replay nonincident pair");
            ++changed;
        }
        if(changed!=1) throw std::runtime_error("replay codimension");
        int found=0;
        for(int a=0;a<3;++a) if(!(fp[a]&1)) for(int s:{-1,1}) {
            auto q=fp;q[a]+=s;
            if(q[a]<0 || q[a]>=g.n[a]) continue;
            Id other=g.id(q[0],q[1],q[2]);
            if(c[other]) {if(other!=p) throw std::runtime_error("replay face not free");++found;}
        }
        if(found!=1) throw std::runtime_error("replay parent count");
        for(int a=0;a<3;++a) if(!(pp[a]&1)) for(int s:{-1,1}) {
            auto q=pp;q[a]+=s;
            if(q[a]>=0 && q[a]<g.n[a] && c[g.id(q[0],q[1],q[2])])
                throw std::runtime_error("replay parent not maximal");
        }
        c[f]=c[p]=0;
    }
}
template<class T> void write_binary(const std::string& path,const std::vector<T>& values) {
    std::ofstream out(path,std::ios::binary);
    if(!out) throw std::runtime_error("cannot open output");
    out.write(reinterpret_cast<const char*>(values.data()),std::streamsize(values.size()*sizeof(T)));
    if(!out) throw std::runtime_error("output write failed");
}
int main(int argc,char** argv) {
    try {
        if(argc!=6) throw std::runtime_error("usage: kernel NX NY NZ MASK PREFIX");
        int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]),nz=std::stoi(argv[3]);
        if(std::min({nx,ny,nz})<1 || std::max({nx,ny,nz})>512) throw std::runtime_error("invalid dimensions");
        Grid g(nx,ny,nz);std::vector<uint8_t> mask(uint64_t(nx)*ny*nz);
        std::ifstream in(argv[4],std::ios::binary);in.read(reinterpret_cast<char*>(mask.data()),mask.size());
        if(in.gcount()!=std::streamsize(mask.size()) || in.peek()!=EOF) throw std::runtime_error("mask size mismatch");
        for(auto x:mask) if(x>1) throw std::runtime_error("nonboolean mask");
        Cells c=build(g,mask,nx,ny,nz);auto initial=counts(g,c);
        std::priority_queue<uint64_t,std::vector<uint64_t>,std::greater<uint64_t>> heap;
        Cells pending(g.size,0);
        auto enqueue=[&](Id f) {
            if(pending[f] || !c[f]) return;
            int64_t p=free_parent(g,c,f);
            if(p>=0) {heap.push(uint64_t(3-g.dim(Id(p)))*g.size+f);pending[f]=1;}
        };
        for(uint64_t f=0;f<g.size;++f) if(c[f] && g.dim(Id(f))<3) enqueue(Id(f));
        std::vector<std::array<Id,2>> moves;
        while(!heap.empty()) {
            Id face=Id(heap.top()%g.size);heap.pop();pending[face]=0;
            int64_t parent=free_parent(g,c,face);if(parent<0) continue;
            Id par=Id(parent);c[face]=c[par]=0;moves.push_back({face,par});
            auto p=g.xyz(par);
            int lo[3],hi[3];for(int a=0;a<3;++a){lo[a]=p[a]-(p[a]&1);hi[a]=p[a]+(p[a]&1);}
            for(int x=lo[0];x<=hi[0];++x) for(int y=lo[1];y<=hi[1];++y) for(int z=lo[2];z<=hi[2];++z) {
                Id f=g.id(x,y,z);if(f!=par) enqueue(f);
            }
        }
        auto terminal_counts=counts(g,c);
        Cells checked=build(g,mask,nx,ny,nz);replay(g,checked,moves);
        if(checked!=c) throw std::runtime_error("replayed endpoint differs");
        std::vector<Id> terminal;for(uint64_t i=0;i<g.size;++i) if(c[i]) terminal.push_back(Id(i));
        std::string prefix=argv[5];write_binary(prefix+".moves",moves);write_binary(prefix+".terminal",terminal);
        std::ofstream out(prefix+".json");
        out<<"{\"replay_valid\":true,\"collapses\":"<<moves.size()<<",\"initial_cell_counts\":[";
        for(int i=0;i<4;++i) out<<(i?",":"")<<initial[i];
        out<<"],\"terminal_cell_counts\":[";for(int i=0;i<4;++i) out<<(i?",":"")<<terminal_counts[i];
        out<<"],\"cell_shape\":["<<g.n[0]<<","<<g.n[1]<<","<<g.n[2]<<"]}\n";
        if(!out) throw std::runtime_error("metadata write failed");
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
