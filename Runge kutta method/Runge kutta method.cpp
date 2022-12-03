//
//  Runge kutta method.cpp
//  Runge kutta method
//
//  Created by WC on 2022/12/3.
//

#include <iostream>
#include <vector>
using namespace std;

template <class state_type>
void resize(const state_type &in, state_type &out) {
    // standard implementation works for containers
    using std::size;
    out.resize(size(in));
}

// specialization for std:: array
template <class T, std:: size_t N>
void resize(const std::array <T, N> &, std::array <T,N>& ) {
    /* arrays don't need resizing */
}

template <class state_type, class time_type>
void component_sum(state_type &w, state_type &u, time_type v) {
    using std::begin, std::end;
    auto w_begin = begin(w); auto w_end = end(w);
    auto u_begin = begin(u);
    for( ; w_begin != w_end ; w_begin++) {
        * w_begin += v * (* u_begin++);
    }
}

template <class state_type, class time_type>
void container_sum(state_type &w, state_type &x, vector<state_type> u, std::vector<time_type> v) {
    using std::begin, std::end;
    state_type w0;
    resize(w, w0);
    auto u_begin = begin(u); auto u_end = end(u);
    auto v_begin = begin(v);
    for(component_sum(w0, x, 1) ; u_begin != u_end ; ) {
        component_sum(w0, *(u_begin++), *(v_begin++));
    }
    w=w0;
}

template <class state_type , class time_type = double>
class runge_kutta4 {
public:
    template <typename System>
    void do_step(System &system , state_type &x, time_type t, time_type dt) {
        adjust_size(x);
        const time_type dt2 = dt/2, dt3 = dt/3, dt6 = dt/6;
        typedef std::vector<state_type> ct_type;
        typedef std::vector<time_type> co_type;

        system(x, s[0] , t);

        container_sum(x_tmp , x , ct_type{s[0]} , co_type{dt2});
        system(x_tmp , s[1] , t + dt2);
        
        container_sum(x_tmp , x , ct_type{s[1]} , co_type{dt2});
        system(x_tmp , s[2] , t + dt2);

        container_sum(x_tmp , x , ct_type{s[2]} , co_type{dt});
        system(x_tmp , s[3] , t + dt);

        container_sum(x , x , s , co_type{dt6, dt3, dt3, dt6});
    }
private:
    state_type x_tmp;
    vector<state_type> s{x_tmp, x_tmp, x_tmp, x_tmp};
    void adjust_size(const state_type &x) {
        resize(x, x_tmp);
        for(int i = 0 ; i < std::size(s) ; i++) {
            resize(x, s[i]);
        }
    }
};

template <class state_type , class time_type = double>
class runge_kutta3 {
public:
    template <typename System>
    void do_step(System &system , state_type &x, time_type t, time_type dt) {
        adjust_size(x);
        const time_type dt2 = dt/2, dt6 = dt/6, dt_=-dt, two_dt = 2*dt, three_dt = 3*dt/2;
        typedef std::vector<state_type> ct_type;
        typedef std::vector<time_type> co_type;

        system(x, s[0] , t);

        container_sum(x_tmp , x , ct_type{s[0]} , co_type{dt2});
        system(x_tmp , s[1] , t + dt2);
        
        container_sum(x_tmp , x , ct_type{s[0], s[1]} , co_type{dt_, two_dt});
        system(x_tmp , s[2] , t + dt);

        container_sum(x , x , s , co_type{dt6, three_dt, dt6});
    }
private:
    state_type x_tmp;
    vector<state_type> s{x_tmp, x_tmp, x_tmp};
    void adjust_size(const state_type &x) {
        resize(x, x_tmp);
        for(int i = 0 ; i < std::size(s) ; i++) {
            resize(x, s[i]);
        }
    }
};
