@extends('layouts.app')

@section('content')
    <div class="container">
        <div class="row justify-content-center">
            <div class="col-md-8">
                <div class="card">
                    <div class="card-header">Change your password</div>
                    <div class="card-body">
                        <form method="post" action="{{ route('do-change-password') }}">
                            @csrf
                            <div class="form-group">
                                <label for="current_password">Current password</label>
                                <input type="password"
                                       class="form-control @error('current_password') is-invalid @enderror"
                                       id="current_password" name="current_password">
                                @error('current_password')
                                <div class="invalid-feedback">{{ $message }}</div>
                                @enderror
                            </div>
                            <div class="form-group">
                                <label for="new_password">New password</label>
                                <input type="password"
                                       class="form-control @error('new_password') is-invalid @enderror"
                                       id="new_password" name="new_password">
                                @error('new_password')
                                <div class="invalid-feedback">{{ $message }}</div>
                                @enderror
                            </div>
                            <div class="form-group">
                                <label for="new_password_confirmation">Confirm password</label>
                                <input type="password"
                                       class="form-control @error('new_password_confirmation') is-invalid @enderror"
                                       id="new_password_confirmation" name="new_password_confirmation">
                                @error('new_password_confirmation')
                                <div class="invalid-feedback">{{ $message }}</div>
                                @enderror
                            </div>

                            <button type="submit" class="btn btn-primary"><i class="fas fa-key"></i>&nbsp;Change Password</button>
                            <a href="{{ route('home') }}" class="btn btn-secondary"><i class="fas fa-arrow-left"></i>&nbsp;Go Back</a>
                        </form>
                    </div>
                </div>
            </div>
        </div>
    </div>
@endsection
