@extends('layouts.app')

@section('content')
    <div class="container">
        <div class="row justify-content-center">
            <div class="col-md-10">
                <div class="card">
                    <div class="card-header">Users</div>
                    <div class="card-body">
                        @if (session('status'))
                            <div class="alert alert-success alert-dismissible fade show" role="alert">
                                {{ session('status') }}
                                <button type="button" class="close" data-dismiss="alert" aria-label="Close">
                                    <span aria-hidden="true">&times;</span>
                                </button>
                            </div>
                        @endif
                        <div class="row">
                            <div class="col-2 offset-10">
                                <a href="{{ route('users-create') }}" class="btn btn-success">
                                    <i class="far fa-plus-square"></i>&nbsp;New user</a>
                            </div>
                        </div>
                        <div class="table-responsive pt-2">
                            <table class="table table-striped">
                                <thead>
                                <tr>
                                    <th scope="col">#</th>
                                    <th scope="col">Name</th>
                                    <th scope="col">E-mail</th>
                                    <th scope="col" class="text-center">Admin</th>
                                    <th scope="col" class="text-center">Actions</th>
                                </tr>
                                </thead>
                                <tbody>
                                @forelse($users as $user)
                                    <tr>
                                        <th scope="row">{{ $user->id }}</th>
                                        <td>{{ $user->name }}</td>
                                        <td>{{ $user->email }}</td>
                                        <td class="text-center">
                                            @if($user->admin)
                                                <i class="fas fa-check-square"></i>
                                            @else
                                                <i class="fas fa-square"></i>
                                            @endif
                                        </td>
                                        <td class="text-center">
                                            <a href="{{ route('users-update', $user) }}" class="btn btn-sm btn-primary">
                                                <i class="fas fa-edit"></i>&nbsp;Edit
                                            </a>
                                            <a href="{{ route('users-delete', $user) }}" class="btn btn-sm btn-danger">
                                                <i class="fas fa-trash-alt"></i>&nbsp;Delete
                                            </a>
                                        </td>
                                    </tr>
                                @empty
                                    <tr>
                                        <td rowspan="5">No users found!</td>
                                    </tr>
                                @endforelse
                                </tbody>
                            </table>
                        </div>
                        <div class="row">
                            <div class="col-10">
                                {{ $users->links() }}
                            </div>
                            <div class="col-2">
                                <a href="{{ route('home') }}" class="btn btn-secondary"><i
                                        class="fas fa-arrow-left"></i>&nbsp;Go Back</a>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
@endsection
